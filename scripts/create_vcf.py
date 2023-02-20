"From a GFA file, tries to create a VCF"
from argparse import ArgumentParser, SUPPRESS
from BubbleGun.Graph import Graph as bGraph
from BubbleGun.find_bubbles import find_bubbles
from gfagraphs import Graph as gGraph
from importlib.metadata import version


def get_all_bubbles(gfa_file: str) -> dict:
    """Gets all bubbles inside a graph

    Args:
        gfa_file (str): a pangenome GFA file

    Returns:
        dict: a description of the bubbles
    """
    find_bubbles(bg_graph := bGraph(gfa_file))
    return {i: [
        val.list_bubble(),
        val.is_insertion(),
        val.is_simple(),
    ] for i, val in enumerate(bg_graph.bubbles.values())}


def get_graph_structure(gfa_file: str, gfa_version: str, reference_path: str, chromosome: str | int) -> dict:
    """Given a GFA graph, extracts variants

    Args:
        gfa_file (str): _description_
        gfa_version (str): _description_
        reference_path (str): _description_
    """
    variants: dict = dict()
    segments: dict = {seg.datas['name']: seg for seg in gGraph(
        gfa_file, gfa_version, with_sequence=True).segments}
    bubbles: dict = get_all_bubbles(gfa_file)
    for id, (node_list, is_insert, is_simple) in bubbles.items():
        if is_insert or is_simple:
            # Simple bubble
            if is_insert:
                # verifies if its the reference that goes through
                _, source, insert = node_list
                if reference_path in (ins := segments[insert]).datas['PO'].keys():
                    # deletion from reference
                    var, ref, path_offset = ins.datas['seq'][
                        0], ins.datas['seq'], ins.datas['PO']
                else:
                    # insertion from reference
                    ref, var, path_offset, ro = ins.datas['seq'][
                        0], ins.datas['seq'], ins.datas['PO'],  segments[source].datas['PO'][reference_path][1]
            else:
                _, _, node_a, node_b = node_list
                low = segments[node_b]
                if reference_path in (high := segments[node_a]).datas['PO'].keys():
                    var, ref, path_offset = low.datas['seq'], high.datas['seq'], {
                        **low.datas['PO'], **high.datas['PO']}
                else:
                    ref, var, path_offset = low.datas['seq'], high.datas['seq'], {
                        **low.datas['PO'], **high.datas['PO']}
            try:
                ref_offset = ro
            except UnboundLocalError:
                ref_offset = path_offset[reference_path][0]
            variants[id] = {'ref': ref, 'alt': var, 'ref_offset': ref_offset, 'alt_offset': [
                (name, alt) for name, alt in path_offset.items() if name != reference_path], "chrom": chromosome, "superbubble": False}
        else:
            # Superbubble, woosh
            for node in node_list[2:]:
                if not reference_path in (current := segments[node]).datas['PO'].keys():
                    # If in reference, not a variant, so exclunding those
                    ref, var, path_offset = current.datas['seq'][0], current.datas['seq'], {
                        **current.datas['PO'], reference_path: segments[node_list[1]].datas['PO'][reference_path]}
                    variants[id] = {'ref': ref, 'alt': var, 'ref_offset': path_offset[reference_path], 'alt_offset': [
                        (name, alt) for name, alt in path_offset.items() if name != reference_path], "chrom": chromosome, "superbubble": True}
    return variants


def render_vcf(output_file: str, graph_variants: dict) -> None:
    "Creates a VCF file"
    vcf_header: str = f"##fileformat=VCFv4.2\n##source=pangraphs v{version('pangraphs')}\n##INFO=<ID=NA,Number=1,Type=String,Description=\"Name of the alternative path\">\n##INFO=<ID=AO,Number=1,Type=Integer,Description=\"Offset (start) for alternative sequence\">\n##INFO=<ID=OR,Number=1,Type=String,Description=\"Orientation for reference sequence\">\n##INFO=<ID=OA,Number=1,Type=String,Description=\"Orientation for alternate sequence\">\n##INFO=<ID=BS,Number=0,Type=Flag,Description=\"If variant is inside a superbubble\">\n#CHROM\tID\tPOS\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    quality: float = 60.0
    filter_info: str = '.'
    with open(output_file, 'w', encoding='utf-8') as vcf_writer:
        vcf_writer.write(vcf_header)
        for id, variant in graph_variants.items():
            for subvar, suboff in variant['alt_offset']:
                flags = ";BS" if variant['superbubble'] else ''
                vcf_writer.write(
                    f"{variant['chrom']}\t{id}\t{variant['ref_offset']}\t{variant['ref']}\t{variant['alt']}\t{quality}\t{filter_info}\tNA={subvar};AO={suboff[0]};OR={variant['ref_offset']};OA={suboff[2]}{flags}\n")


if __name__ == "__main__":

    parser = ArgumentParser(add_help=False)
    parser.add_argument("file", type=str, help="Input GFA file")
    parser.add_argument("output", type=str,
                        help="Output path for VCF (with extension)")
    parser.add_argument(
        "-g", "--gfa_version", help="Tells the GFA input style", required=True, choices=['rGFA', 'GFA1', 'GFA1.1', 'GFA1.2', 'GFA2'])
    parser.add_argument('-h', '--help', action='help', default=SUPPRESS,
                        help='This script aims to extract variants from a GFA graph.')
    parser.add_argument(
        "-c", "--chromosom", help="Name of the chromosom the graph is from", required=True, type=str)
    parser.add_argument(
        "-r", "--reference", help="Name of the reference path inside graph", required=True, type=str)
    args = parser.parse_args()

    render_vcf(args.output, get_graph_structure(
        args.file, args.gfa_version, args.reference, args.chromosom))
