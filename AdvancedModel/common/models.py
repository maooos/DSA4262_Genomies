"""Three permutation-invariant read-bag classifiers with site-level supervision."""

from dataclasses import asdict, dataclass
import math

import torch
from torch import nn
from torch.nn import functional as F


@dataclass
class ModelConfig:
    architecture: str = "noisy_or"
    embedding_dim: int = 8
    hidden_dim: int = 128
    read_dim: int = 64
    attention_dim: int = 64
    dropout: float = 0.1
    heads: int = 4
    set_layers: int = 2

    def __post_init__(self):
        if self.architecture not in {"noisy_or", "gated_attention", "set_transformer"}:
            raise ValueError("Unknown MIL architecture.")
        if min(self.embedding_dim, self.hidden_dim, self.read_dim, self.attention_dim, self.heads, self.set_layers) < 1:
            raise ValueError("Model dimensions must be positive.")
        if not 0 <= self.dropout < 1:
            raise ValueError("dropout must lie in [0, 1).")
        if self.architecture == "set_transformer" and self.read_dim % self.heads:
            raise ValueError("read_dim must be divisible by heads.")


class GatedAttention(nn.Module):
    def __init__(self, dimension, attention_dimension):
        super().__init__()
        self.value = nn.Linear(dimension, attention_dimension)
        self.gate = nn.Linear(dimension, attention_dimension)
        self.score = nn.Linear(attention_dimension, 1, bias=False)

    def forward(self, reads, mask):
        logits = self.score(torch.tanh(self.value(reads)) * torch.sigmoid(self.gate(reads))).squeeze(-1)
        weights = logits.masked_fill(~mask, -torch.inf).softmax(dim=1)
        return torch.sum(weights.unsqueeze(-1) * reads, dim=1), weights


class SetAttentionBlock(nn.Module):
    """Self-attention block (SAB); no ordering/position embedding among reads."""
    def __init__(self, dimension, heads, dropout):
        super().__init__()
        self.norm1, self.norm2 = nn.LayerNorm(dimension), nn.LayerNorm(dimension)
        self.attention = nn.MultiheadAttention(dimension, heads, dropout=dropout, batch_first=True)
        self.feedforward = nn.Sequential(nn.Linear(dimension, dimension * 2), nn.GELU(),
                                        nn.Dropout(dropout), nn.Linear(dimension * 2, dimension))
        self.dropout = nn.Dropout(dropout)

    def forward(self, reads, mask):
        normalized = self.norm1(reads)
        attended, _ = self.attention(normalized, normalized, normalized,
                                     key_padding_mask=~mask, need_weights=False)
        reads = reads + self.dropout(attended)
        reads = reads + self.dropout(self.feedforward(self.norm2(reads)))
        return reads.masked_fill(~mask.unsqueeze(-1), 0.)


class MILModel(nn.Module):
    def __init__(self, config):
        super().__init__()
        self.config = config
        self.kmer_embedding = nn.Embedding(1024, config.embedding_dim)
        self.read_encoder = nn.Sequential(
            nn.Linear(15 + 3 * config.embedding_dim, config.hidden_dim),
            nn.LayerNorm(config.hidden_dim), nn.GELU(), nn.Dropout(config.dropout),
            nn.Linear(config.hidden_dim, config.read_dim), nn.LayerNorm(config.read_dim), nn.GELU(),
        )
        self.set_encoder = nn.ModuleList([
            SetAttentionBlock(config.read_dim, config.heads, config.dropout)
            for _ in range(config.set_layers if config.architecture == "set_transformer" else 0)
        ])
        if config.architecture == "noisy_or":
            self.read_head = nn.Linear(config.read_dim, 1)
            # A low initial read probability avoids immediate Noisy-OR saturation.
            nn.init.constant_(self.read_head.bias, math.log(0.002 / 0.998))
        else:
            self.pooling = GatedAttention(config.read_dim, config.attention_dim)
            self.site_head = nn.Sequential(nn.Linear(config.read_dim, 32), nn.GELU(),
                                           nn.Dropout(config.dropout), nn.Linear(32, 1))
        self.fraction_head = nn.Linear(config.read_dim, 1)

    def forward(self, features, tokens, mask):
        if features.ndim != 3 or features.shape[-1] != 15 or mask.shape != features.shape[:2]:
            raise ValueError("Expected [batch, reads, 15] features and a matching mask.")
        if not mask.any(dim=1).all():
            raise ValueError("Every site must have at least one real read.")
        context = self.kmer_embedding(tokens).flatten(start_dim=1)
        context = context.unsqueeze(1).expand(-1, features.shape[1], -1)
        clean_features = features.masked_fill(~mask.unsqueeze(-1), 0.)
        reads = self.read_encoder(torch.cat([clean_features, context], dim=-1))
        reads = reads.masked_fill(~mask.unsqueeze(-1), 0.)
        for block in self.set_encoder:
            reads = block(reads, mask)
        if self.config.architecture == "noisy_or":
            read_logits = self.read_head(reads).squeeze(-1)
            # log(1-P_site) = sum log(1-p_read); stable even for many reads.
            log_survival = (-F.softplus(read_logits) * mask).sum(dim=1)
            log_modified = torch.log((-torch.expm1(log_survival)).clamp_min(1e-12))
            logits = log_modified - log_survival
            weights = mask.to(reads.dtype) / mask.sum(dim=1, keepdim=True)
            pooled = torch.sum(reads * weights.unsqueeze(-1), dim=1)
        else:
            pooled, weights = self.pooling(reads, mask)
            logits = self.site_head(pooled).squeeze(-1)
        return {"logits": logits, "fraction_logits": self.fraction_head(pooled).squeeze(-1),
                "attention": weights}

    def specification(self):
        return {**asdict(self.config), "parameters": sum(p.numel() for p in self.parameters())}


def site_loss(outputs, batch, dataset_weights=(1., 1., 1.), auxiliary_weight=0., positive_weight=1.):
    """One supervised loss per site; no propagation of site labels to reads."""
    logits = outputs["logits"]
    weights = logits.new_tensor(dataset_weights)[batch["dataset"]]
    binary = F.binary_cross_entropy_with_logits(
        logits, batch["label"], reduction="none", pos_weight=logits.new_tensor(positive_weight),
    )
    loss = (binary * weights).sum() / weights.sum()
    synthetic = batch["dataset"] == 2
    if auxiliary_weight and synthetic.any():
        # Mixture fraction is a site-level auxiliary target, not a read label.
        auxiliary = (torch.sigmoid(outputs["fraction_logits"][synthetic]) - batch["fraction"][synthetic]).square()
        loss = loss + auxiliary_weight * (auxiliary * weights[synthetic]).sum() / weights.sum()
    return loss
