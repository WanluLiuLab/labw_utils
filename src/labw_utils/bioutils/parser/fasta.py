def extract_fasta_name(line: str, full_header: bool) -> str:
    if full_header:
        return line[1:].strip()
    else:
        return line[1:].strip().split(' ')[0].split('\t')[0]
