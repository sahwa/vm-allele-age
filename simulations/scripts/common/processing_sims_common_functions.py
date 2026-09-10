def write_plink_bed(dosages, prefix, positions, ref, alt, iids, chrom=1):
    """
    dosages : (n_var, n_ind) uint8, count of the ALT (derived) allele, 0/1/2
    prefix  : output path prefix

    A1 is set to ALT and A2 to REF, matching what `plink2 --vcf` does by
    default, so allele orientation is unchanged from the old VCF route.

    PLINK1 .bed: 3 magic bytes, then 2 bits per genotype, 4 per byte, with the
    first individual in the lowest-order bits.
        00 = hom A1, 01 = missing, 10 = het, 11 = hom A2
    With A1 = ALT: dosage 2 -> 00 (0), dosage 1 -> 10 (2), dosage 0 -> 11 (3)
    """
    n_var, n_ind = dosages.shape
    n_bytes = (n_ind + 3) // 4
    pad = n_bytes * 4 - n_ind

    lut = np.array([3, 2, 0], dtype=np.uint8)   # indexed by dosage

    with open(f"{prefix}.bed", "wb") as f:
        f.write(bytes([0x6C, 0x1B, 0x01]))      # magic + SNP-major
        for start in range(0, n_var, 2000):
            block = lut[dosages[start:start + 2000]]
            if pad:
                block = np.pad(block, ((0, 0), (0, pad)))
            block = block.reshape(block.shape[0], n_bytes, 4)
            packed = (block[:, :, 0]
                      | (block[:, :, 1] << 2)
                      | (block[:, :, 2] << 4)
                      | (block[:, :, 3] << 6)).astype(np.uint8)
            packed.tofile(f)

    # IDs match the old `--set-all-var-ids '@:#:$r:$a'` scheme
    snp_ids = [f"{chrom}:{int(p)}:{r}:{a}" for p, r, a in zip(positions, ref, alt)]
    pd.DataFrame({
        "chr": chrom,
        "snpid": snp_ids,
        "cm": 0,
        "pos": positions.astype(np.int64),
        "a1": alt,
        "a2": ref,
    }).to_csv(f"{prefix}.bim", sep="\t", index=False, header=False)

    pd.DataFrame({
        "fid": 0, "iid": iids, "pid": 0, "mid": 0, "sex": 0, "pheno": -9
    }).to_csv(f"{prefix}.fam", sep="\t", index=False, header=False)


