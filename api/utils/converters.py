
def convert_variant_id(variant_id):
    chr, pos, allele_part = variant_id.split("_")
    ref, alt = allele_part.split("/")
    pos = int(pos)
    return chr, pos, ref, alt