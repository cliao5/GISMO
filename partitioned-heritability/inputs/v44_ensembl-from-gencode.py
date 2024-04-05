
with open("gencode.v44.chr_patch_hapl_scaff.basic.annotation.gff3", 'r') as r:
    with open("gencode.v44.chr_patch_hapl_scaff.basic.annotation.gff3.genes", 'w') as w:
        w.write("seqnames\tstart\tend\tgene_id\tgene_name\n")
        # Iterate through .gff3 file line by line
        while True:
            line = r.readline()
            if not line:
                break

            # Skip header
            if line[0] == "#":
                continue

            annotations = line.split()

            # Skip annotations that aren't "gene" feature types
            if annotations[2] != "gene":
                continue

            # Otherwise, collect details that we want
            chr = annotations[0].replace("chr", "")
            start = annotations[3]
            end = annotations[4]
            gene_name = ""
            gene_id = ""

            for a in annotations[8].split(";"):
                b = a.split("=")
                
                if b[0] == "gene_name":
                    gene_name = b[1]
                
                elif b[0] == "gene_id":
                    gene_id = b[1].split(".")[0]
            
            # Write only if we have everything
            if chr and start and end and gene_name and gene_id:
                w.write(f"{chr}\t{start}\t{end}\t{gene_id}\t{gene_name}\n")



        