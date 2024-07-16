# Commands description
HELP_COMMAND_ISOLATE: str = "Isolates a subgraph within a graph in-between two positions, according to a path or walk in the graph. Extraction only considers other paths and walks that cross both positions."
HELP_COMMAND_OFFSET: str = "Add path offsets to a graph. Adds a JSON string, PO (Path Offset) positions, relative to paths in the GFA file. Hence, PO:J:{'w1':[(334,336,'+')],'w2':[(245,247,'-'),(336,338,'-')]} tells that the walk/path w1 contains the sequence starting at position 334 and ending at position 336, and the walk/path w2 contains the sequence starting at the offset 245 (ending 247) and crosses it a second time between position 336 and 338, and that the sequences are reversed one to each other. Note that any non-referenced walk in this field means that the node is not inside the given walk."
HELP_COMMAND_COMPLETE: str = "Asserts if the graph is a complete pangenome graph. A complete pangenome graph is a sequence graph where all the genomes used to build the graph are exactly contained in the graph."
HELP_COMMAND_GRAPHER: str = "Creates a html view of the graph. Huge graphs may take long time to display, or might be messy. Advice would be to extract parts you want to display (with pancat isolate for instance), then computes the vizualisation on the selected part."
HELP_COMMAND_MULTIGRAPHER: str = "Creates a html view of the alignment of two graphs. Huge graphs may take long time to display, or might be messy. Advice would be to extract parts you want to display (with pancat isolate for instance), then computes the vizualisation on the selected part."
HELP_COMMAND_STATS: str = "Retrieves basic stats on a pangenome graph."
HELP_COMMAND_RECONSTRUCT: str = "Reconstruct linear sequences from a GFA graph. Given a GFA file with paths, reconstruct a linear sequence for each haplotype between two offsets."
HELP_COMMAND_UNFOLD: str = "[WIP] Break cycles of a graph with paths or walks. Aims to linearize the graph by creating new nodes and connecting them in the graph."
HELP_COMMAND_COMPRESS: str = "[WIP] Does a compression of the graph, merging simple non-conflicting substitution bubbles."
HELP_COMMAND_EDIT: str = "Aligns graphs and tires to compute the minimal set of events that explains how we go from one graph to another."

# Inputs
HELP_INPUT_FILE_WITH_PATHS: str = "Path to a GFA-like file. GFA must contain P or W lines."
HELP_INPUT_EDITION_FILE: str = "Path-level .json edition file created with pancat edit between file_A and file_B"
HELP_INPUT_FILE_GENERAL: str = "Path to a GFA or rGFA-like file."
HELP_INPUT_REFERENCE: str = "The path used for positions. Extraction happen on this sequence, in-between the two specified positions. The two nodes the positions are included into are analysed (resp. named source and sink), and a intersection of sets of the haplotypes going through both nodes is used for extraction, meaning the only paths that are extracted along the reference are the ones that goes trough both source and sink nodes. For W-lines, it is the 4th field, whereas for P-lines it is the 2nd one."
HELP_INPUT_START_POSITION: str = "The starting point on --reference where the extraction should start, in basepairs from the start of the selected path."
HELP_INPUT_END_POSITION: str = "The ending point on --reference where the extraction should end, in basepairs from the start of the selected path."
HELP_INPUT_PIPELINE: str = "Tab-separated mapping between path names and path to files."

# Outputs
HELP_OUTPUT_FILE_GFA: str = "Output path for the graph. It should be provided with the .gfa extension."
HELP_OUTPUT_FILE_HTML: str = "Output path for the visualisation. It should be provided with the .html extension."
HELP_OUTPUT_FOLDER: str = "Path to a output folder. The path could be given with trailing bacslash but should'nt contain any extension."
HELP_OUTPUT_JSON_LOG: str = "Path to a .json or .log output file for results. If container is a .json, the editions are stored in a dict structure but it requires memory to handle them all at a time. As a .log file, each path gets its own container, and memory is discarded and allocated on-the-fly, which is far better for huge graphs."

# Parameters
HELP_PARAM_BOUNDARIES: str = "One or more integers to use as boundaries for node size classes (e.g. -b 50 2000 will set 3 classes : one for nodes in range 0-50bp, one for nodes in range 51-2000 bp and one for nodes in range 2001-inf bp)."
HELP_PARAM_MINIMIZE: str = "Saves the graph with a minimal set of informations (for compatibiliy purposes with software such as odgi or vg)"
HELP_PARAM_SPLIT: str = "If enlabed, each record will be "
HELP_PARAM_LENGTH: str = "Maximum length of substitutions to compress"
HELP_PARAM_PATHSELECT: str = "Name(s) for the paths you want to apply method onto. For W-lines, it is the 4th field, whereas for P-lines it is the 2nd one."
HELP_PARAM_REGEXP: str = "Regexp to filer on for path names. Dealults to no filtering."
HELP_PARAM_GRAPH_LEVEL: str = "Asks to perform edition computation at graph level."
HELP_PARAM_THREADS: str = "Number of threads used for parallelization."
HELP_PARAM_VERBOSE: str = "Print to log extended information about what is going on."
