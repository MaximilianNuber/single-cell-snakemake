from __future__ import annotations
from typing import Protocol, Any, runtime_checkable
import numpy as np
from .types import Result

@runtime_checkable
class EmbeddingComputer(Protocol):
    def compute(self, X_counts: Any, obs: Any, var: Any) -> Result: ...

@runtime_checkable
class NeighborGraphBuilder(Protocol):
    def build(self, X_embed: np.ndarray) -> Result: ...

@runtime_checkable
class ManifoldEmbedder(Protocol):
    def embed(self, connectivities, distances=None) -> Result: ...