# Recipe: spatial transcriptomics

1. Provide a Seurat spatial RDS or Space Ranger directory.
2. Preserve spatial coordinates, image and `spatial_sample_id`.
3. Normalize and cluster each supported object.
4. Identify spatially variable genes.
5. Provide a compatible annotated scRNA reference for transfer/deconvolution.
6. Run label transfer and optional RCTD only when reference gates pass.
7. Run neighborhood permutation analysis without erasing slice identity.
8. Export one independent map per slice, cell type, score or gene.
