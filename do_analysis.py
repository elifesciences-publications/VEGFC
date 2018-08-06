#!/usr/bin/python3
# -*- coding: UTF-8 -*-
#
#
# UBUNTU PACKAGES
#
# This script has been developed to run under Linux and was specifically tested ONLY under Ubuntu 16.04 and 18.04.
# It requires the installation of the following ubuntu packages:
#
# python3-setuptools
# clustalw
# mafft
# dialign-tx
# poa
# probcos
# muscle
# kalign
# amap-align
# proda
# prank
# t-coffee
# phyml
#
# sudo apt -y --show-progress install inkscape python3-setuptools python3-pyqt4 python3-pyqt4.qtopengl python3-pip autoconf t-coffee clustalw mafft dialign-tx poa probcons muscle kalign amap-align proda prank phyml t-coffee imagemagick build-essential libblas-dev liblapack-dev zlib1g-dev libcairo2-dev libcurl4-openssl-dev python3-numpy python3-lxml python3-six
#
#
# For some reason the Ubuntu package names some of the alignment executables differently and some links need to be
# created in order for t-coffee to find them:
#
# sudo ln -s /usr/bin/dialign-tx /bin/dialign-t
# sudo ln -s /usr/bin/clustalw /bin/clustalw2
#
#
# MANUAL INSTALLATIONS
#
# In addition, the pcma executable need to be manually downloaded from http://prodata.swmed.edu/download/pub/PCMA/pcma.tar.gz
# and compiled because it is not available from the Ubuntu repository:
#
# wget http://prodata.swmed.edu/download/pub/PCMA/pcma.tar.gz
# tar -xvzf pcma.tar.gz
# cd pcma
# make
# sudo cp pcma /bin
#
#
# OPTIONAL INSTALLS
#
# The multicore-enabled version of phyml (phyml-mpi) is not available as a precompiled Ubuntu package and needs
# to be installed manually, but they single-core version works as well (is just slower). The command to execute
# the multicore version is: "mpirun -n 4 phyml-mpi -i " + PHYLIP_ALIGNED_TRIMMED_CODED +  " -d aa -b -1"
# If this script is run on a (headless) server, the xvfb package is required since the ete3 package requires the presence of x.org.
#
#
# PYTHON MODULES
#
# The following Python modules need to be installed:
#
# biopython
# ete3
#
# sudo pip3 install biopython ete3

import argparse, subprocess, Bio, os, sys
#from Bio.Blast import NCBIWWW
#from Bio.Blast import NCBIXML
#from Bio.Blast.Applications import NcbipsiblastCommandline
from Bio import Entrez
from Bio import SeqIO
from ete3 import Tree, TreeStyle, TextFace, NodeStyle, SequenceFace, ImgFace, SVGFace, faces, add_face_to_node

FASTA = 'sequences.fasta'
FASTA_ALIGNED = 'sequences_aligned.fasta'
FASTA_ALIGNED_TRIMMED = 'sequences_aligned_trimmed.fasta'
FASTA_ALIGNED_TRIMMED_CODED = 'sequences_aligned_trimmed_coded.fasta'
PHYLIP_ALIGNED_TRIMMED_CODED = 'sequences_aligned_trimmed_coded.phy'
PHYLIP_ALIGNED_TRIMMED_DECODED = 'sequences_aligned_trimmed_decoded.newick'
EVALUATIONFILE = 'sequences_aligned_trimmed_score.fasta_aln'
IMG_PATH = "./images/"
SVG_TREEFILE = 'sequences_aligned.svg'
REFERENCE_SEQUENCE = "NP_005420.1" # Homo sapiens


parser = argparse.ArgumentParser()
parser.add_argument("options", help = "Possible command line arguments: download - only download sequences; align - only align sequences; evaluate - only evaluate sequence alignment; drawtree - only draw tree", nargs = '*', default = ['download', 'align', 'evaluate', 'drawtree'])
args = parser.parse_args()

# Animal SVG silouette images from http://phylopic.org (public domain) or own creations
# It would be nice to have a lungfish VEGF-C, but the genomes have not been made publicly available!
# Protopterus annectens
# Neoceratodus forsteri
# Lepidosiren paradoxa
#
sequence_dictionary =  {'XP_022098898.1':["Starfish.svg", "Starfish", "Acanthaster planci", "Asteroidea"],
         'XP_007894115.1':["Callorhinchus_milii.svg", "Ghost shark", "Callorhinchus milii", "Chondrichthyes"],
#         'NP_002599.1':["Homo_sapiens.svg", "Human PDGF-B", "Homo sapiens", "Mammalia"],
         'XP_020376152.1':["Rhincodon_typus.svg", "Whale shark", "Rhincodon typus", "Chondrichthyes"],
         'ENSEBUT00000000354.1':["Eptatretus_burgeri.svg", "Hagfish", "Eptatretus burgeri", "Myxini", '''>ENSEBUT00000000354.1 PREDICTED: VEGF-C Inshore hagfish [Eptatretus burgeri]
LAIDVLHLHIHPDYLQDNEDIQTDHDPWEIIDTDTFPKGALGPKRIERLTRRLLAASSVD
DLLTLLYPWPEEATAQRCRRGHRTEPQFQAAVVNINWEAIELEWSNTLCAPRQACVPTGP
DSHSVERSLHYRPPCVSLHRCTGCCNDPRRSCTSTAVQHVSKTVIEISLFPELVIRPVTI
SYKNHTECHCLTIPFHNVRPPRSVSKTWRDGGQCGPTSGSCAKGTSWNVEACRCVAQQGV
GEVCGPGMKWNEEMCNCVCWRVCPRGQRLHTQSCGCECALNTRDCFLRARRFDRRKCRCV
TAPCPGAPEVCPVGLGFSEELCRCVPQDWIQGLQRNGG\n'''],
         'LS-transcriptB2-ctg17881':['Leucoraja_erinacea.svg', 'Skate', "Leucoraja erinacea", "Chondrichthyes", '''>LS-transcriptB2-ctg17881 PREDICTED: VEGF-C Little skate [Leucoraja erinacea]
RDQAHSQGQATSQLEQQLRSAASIIELMDIFYPEYRRIQECLQRRSTMAKHARREVEEEQ
EEEEEEEWTEAAAFTVLWREEDLRNIELEWERTQCKPREVCLDLGRELGTATNNFYKPPC
VSVHRCGGCCNNEGFQCINVSTAFVSKTLMEITIPQVGLSRPVVISFINHTACGCHPRHI
FSHSHSIIRRSFHVSPTSCVMGNETCPRGHHWDPHHCGCVSVHEVAAPPASTAEPDVTEG
EFDDFCGPYMVFDEDSCSCVCTNRPSSCHPSKEFDENTCRCVCFNRQHRGLCREEEQEEW
DDDACQCVCRKSCPRHLPLNTNTCTCECSESPASCFRRGKKFDPYTCRCYRLPC\n'''],
         'XP_006632034.2':["Spotted_gar.svg", "Spotted gar", "Lepisosteus oculatus", "Actinopterygii"],
         'NP_001167218.1':["Salmo_salar.svg", "Atlantic salmon", "Salmo salar", "Actinopterygii"],
         'NP_991297.1':["Danio_rerio.svg", "Zebrafish", "Danio rerio", "Actinopterygii"],
         'XP_011610643.1':["Takifugu_rubripes.svg", "Fugu fish", "Takifugu rubripes", "Actinopterygii"],
         'XP_020464669.1':["Monopterus_albus.svg", "Swamp eel", "Monopterus albus", "Actinopterygii"],
         'XP_015809835.1':["Nothobranchius_furzeri.svg", "Killifish", "Nothobranchius furzeri", "Actinopterygii"],
         'XP_023189044.1':["Platyfish.svg", "Platy fish", "Xiphophorus maculatus", "Actinopterygii"],
         'XP_007564695.1':["Poecilia_formosa.svg", "Amazon molly", "Poecilia formosa", "Actinopterygii"],
         'XP_006006690.1':["Coelacant.svg", "Coelacant", "Latimeria chalumnae", "Sarcopterygii"],
         'XP_002933363.1':["Xenopus_tropicalis.svg", "Xenopus", "Xenopus tropicalis", "Amphibia"],
         'XP_018419054.1':["Nanorana_parkeri.svg", "Tibet frog", "Nanorana parkeri", "Amphibia"],
         'XP_015283812.1':["Gekko_japonicus.svg", "Gekko", "Gekko japonicus", "Reptilia"],
         'ETE60014.1':["Ophiophagus_hannah.svg", "King cobra", "Ophiophagus hannah", "Reptilia"],
         'XP_003221689.1':["Lizard.svg", "Anole lizard", "Anolis carolinensis", "Reptilia"],
         'XP_005304228.1':["Turtle.svg", "Painted turtle", "Chrysemys picta bellii", "Reptilia"],
         'XP_006276984.1':["Alligator_mississippiensis.svg", "American Alligator", "Alligator mississippiensis", "Reptilia"],
         'XP_009329004.1':["Pygoscelis_adeliae.svg", "Penguin", "Pygoscelis adeliae", "Aves"],
         'XP_013045797.1':["Anser_cygnoides_domesticus.svg", "Goose", "Anser cygnoides domesticus", "Aves"],
         'XP_420532.3':["Gallus_gallus.svg", "Chicken", "Gallus gallus", "Aves"],
         'XP_009486688.1':["Pelecanus_crispus.svg", "Pelican", "Pelecanus crispus", "Aves"],
         'XP_008490178.1':["Calypte_anna.svg", "Hummingbird", "Calypte anna", "Aves"],
         'XP_009564005.1':["Cuculus_canorus.svg", "Cuckoo", "Cuculus canorus", "Aves"],
         'XP_002189592.1':["Taeniopygia_guttata.svg", "Zebra finch", "Taeniopygia guttata", "Aves"],
         'XP_003415871.1':["Loxodonta_africana.svg", "African elephant", "Loxodonta africana", "Mammalia"],
         'XP_540047.2':["Canis_familiaris.svg", "Dog", "Canis lupus familiaris", "Mammalia"],
         'XP_526740.1':["Pan_troglodytes.svg", "Chimpanzee", "Pan troglodytes", "Mammalia"],
         'NP_005420.1':["Homo_sapiens.svg", "Human", "Homo_sapiens", "Mammalia"],
         'NP_776913.1':["Bos_taurus.svg", "Cattle", "Bos taurus", "Mammalia"],
         'XP_019777186.1':["Tursiops_truncatus.svg", "Dolphin", "Tursiops truncatus", "Mammalia"],
         'XP_004280970.1':["Orcinus_orca.svg", "Orca", "Orcinus orca", "Mammalia"],
         'NP_033532.1':["Mus_musculus.svg", "Mouse", "Mus musculus", "Mammalia"],
         'NP_446105.1':["Rattus_norvegicus.svg", "Rat", "Rattus norvegicus", "Mammalia"],
         'XP_007496150.2':["Monodelphis_domestica.svg", "Opossum", "Monodelphis domestica", "Mammalia"],
         'XP_004465018.1':["Dasypus_novemcinctus.svg", "Armadillo", "Dasypus novemcinctus", "Mammalia"],
         'XP_002709527.1':["Oryctolagus_cuniculus.svg", "Rabbit", "Oryctolagus cuniculus", "Mammalia"],
         'XP_017897598.1':["Capra_hircus.svg", "Goat", "Capra hircus", "Mammalia"]}

def print_subprocess_result(name, out, err):
    if out.decode('utf-8') != '':
        print("Output from " + name + ": " + str(out))
#    if err != None:
    if err.decode('UTF-8') != '':
        print("Error from " + name + ": " + str(err))

def execute_subprocess(comment, bash_command):
    print("\n" + comment, bash_command)
    process = subprocess.Popen(bash_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    output, error = process.communicate()
    process_status = process.wait()
    if output.decode('utf-8') != '':
        print("Output: " + str(output))
    if error.decode('UTF-8') != '':
        print("Error: " + str(error))

def download():
    # Get all fasta sequences from the SEQUENCE_LIST via Entrez
    Entrez.email = "michael@jeltsch.org"
    Entrez.tool = "local_script_under_development"
    file = open(FASTA,"w")
    print("Retrieving sequences from Entrez:\n")
    for item in sequence_dictionary:
        # If the sequence is given locally
        if len(sequence_dictionary[item]) > 4:
            print(sequence_dictionary[item][4])
            file.write(sequence_dictionary[item][4])
        else:
            print("Retrieving " + item)
            try:
                with Entrez.efetch(db="protein", rettype="fasta", retmode="text", id=item) as handle:
                    seq_record = SeqIO.read(handle, "fasta")
                    file.write(seq_record.format("fasta"))
            except Exception as ex:
                print("Problem contacting Blast server. Skipping " + item + ". Error: " + str(ex))
    file.close()

def align():
    execute_subprocess(
        "Generating multiple sequence alignment with the following command:",
        "t_coffee " + FASTA + " -outfile " + FASTA_ALIGNED + " -output=fasta_aln -mode mcoffee")

    execute_subprocess(
        "Trimming multiple sequence alignment with the following command:",
        "t_coffee -other_pg seq_reformat -in " + FASTA_ALIGNED + " -action +extract_block " + REFERENCE_SEQUENCE + " 55 227 > " + FASTA_ALIGNED_TRIMMED)

    execute_subprocess(
        "Converting fasta descriptions part 1 (creating code list) with t_cofeee using the following command:",
        "t_coffee -other_pg seq_reformat -in " + FASTA_ALIGNED_TRIMMED + " -output code_name > code_names.list")

    execute_subprocess(
        "Converting fasta descriptions part 2 (replacing fasta descriptions with codes) with t_cofeee using the following command:",
        "t_coffee -other_pg seq_reformat -code code_names.list -in " + FASTA_ALIGNED_TRIMMED + " > " + FASTA_ALIGNED_TRIMMED_CODED)

    execute_subprocess(
        "Convert into phylip using the following command:",
        "t_coffee -other_pg seq_reformat -in " + FASTA_ALIGNED_TRIMMED_CODED + " -output phylip_aln > " + PHYLIP_ALIGNED_TRIMMED_CODED)

    execute_subprocess(
        "Make tree with the following command:",
        "phyml -i " + PHYLIP_ALIGNED_TRIMMED_CODED +  " -d aa -b -1")

    # phyml adds or doesn't add the .txt extension to the output file (depending on the version) and we need to check for this!
    phyml_output_file = PHYLIP_ALIGNED_TRIMMED_CODED + "_phyml_tree"
    if os.path.isfile(phyml_output_file):
        os.rename(phyml_output_file, phyml_output_file + ".txt")
    execute_subprocess(
        "Decoding tree file file into human-readable format using the following command:",
        "t_coffee -other_pg seq_reformat -decode code_names.list -in " + phyml_output_file + ".txt > " + PHYLIP_ALIGNED_TRIMMED_DECODED)

def evaluate():
    execute_subprocess(
        "Evaluating existing MSA using the following command:",
        "t_coffee -infile " + FASTA_ALIGNED_TRIMMED + " -evaluate -output score_ascii")

    execute_subprocess(
        "Converting ascii score into (pseudo)fasta score using the following command:",
        "t_coffee -other_pg seq_reformat -in sequences_aligned_trimmed.score_ascii -input number_aln -output fasta_aln > " + EVALUATIONFILE)

def drawtree():
    # LOAD FASTA SEQUENCES
    all_sequences = list(SeqIO.parse(FASTA_ALIGNED_TRIMMED, "fasta"))
    dict_of_sequences = {}
    for seq_record in all_sequences:
        dict_of_sequences[seq_record.id] = seq_record.seq
        print(seq_record.id + " : " + seq_record.seq)

    # LOAD TCS EVALUATION FILE
    all_evaluations = list(SeqIO.parse(EVALUATIONFILE, "fasta"))
    dict_of_evaluations = {}
    for seq_record in all_evaluations:
        dict_of_evaluations[seq_record.id] = seq_record.seq
        print(seq_record.id + " : " + seq_record.seq)

    # AMINO ACID COLORING
    amino_acid_fgcolor_dict = {   'A': 'Black',
                                'C': 'Black',
                                'D': 'Black',
                                'E': 'Black',
                                'F': 'Black',
                                'G': 'Black',
                                'H': 'Black',
                                'I': 'Black',
                                'K': 'Black',
                                'L': 'Black',
                                'M': 'Black',
                                'N': 'Black',
                                'P': 'Black',
                                'Q': 'Black',
                                'R': 'Black',
                                'S': 'Black',
                                'T': 'Black',
                                'V': 'Black',
                                'W': 'Black',
                                'Y': 'Black',
                                '-': 'Black'}

    # AMINO ACID BACKGROUND COLORING
    amino_acid_bgcolor_dict = {   'A': 'White',
                                'C': 'White',
                                'D': 'White',
                                'E': 'White',
                                'F': 'White',
                                'G': 'White',
                                'H': 'White',
                                'I': 'White',
                                'K': 'White',
                                'L': 'White',
                                'M': 'White',
                                'N': 'White',
                                'P': 'White',
                                'Q': 'White',
                                'R': 'White',
                                'S': 'White',
                                'T': 'White',
                                'V': 'White',
                                'W': 'White',
                                'Y': 'White',
                                '-': 'White'}

    outgroup_name = 'XP_022098898.1'

    t = Tree(PHYLIP_ALIGNED_TRIMMED_DECODED)
    ts = TreeStyle()
    ts.show_leaf_name = False
    # Zoom in x-axis direction
    ts.scale = 40
    # This makes all branches the same length!!!!!!!
    ts.force_topology = True
    #Tree.render(t, "final_tree_decoded.svg")
    t.set_outgroup(t & outgroup_name)
    ts.show_branch_support = True
    ts.show_branch_length = True
    ts.draw_guiding_lines = True
    ts.branch_vertical_margin = 10 # 10 pixels between adjacent branches
    print(t)

    # Define node styles for different animal classes
    # Mammals
    mammalia = NodeStyle()
    mammalia["bgcolor"] = "Chocolate"
    # Reptiles
    reptilia = NodeStyle()
    reptilia["bgcolor"] = "Olive"
    # Cartilaginous fish
    chondrichthyes = NodeStyle()
    chondrichthyes["bgcolor"] = "SteelBlue"
    # Ray-fenned fish
    actinopterygii = NodeStyle()
    actinopterygii["bgcolor"] = "CornflowerBlue"
    # Lobe-finned fish
    sarcopterygii = NodeStyle()
    sarcopterygii["bgcolor"] = "DarkCyan"
    # Birds
    aves = NodeStyle()
    aves["bgcolor"] = "DarkSalmon"
    # Amphibia
    amphibia = NodeStyle()
    amphibia["bgcolor"] = "DarkSeaGreen"
    # Myxini
    myxini = NodeStyle()
    myxini["bgcolor"] = "LightBlue"

    general_leaf_style = NodeStyle()
    # size of the blue ball
    general_leaf_style["size"] = 15

    # Draws nodes as small red spheres of diameter equal to 10 pixels
    nstyle = NodeStyle()
    nstyle["shape"] = "sphere"
    nstyle["size"] = 15
    nstyle["fgcolor"] = "darkred"
    # Gray dashed branch lines
    #nstyle["hz_line_type"] = 1
    #nstyle["hz_line_color"] = "#cccccc"

    # Applies the same static style to all nodes in the tree if they are not leaves.
    # Note that if "nstyle" is modified, changes will affect to all nodes
    # Apply a separate style to all leaf nodes
    for node in t.traverse():
        if node.is_leaf():
            print("Setting leaf style for node " + node.name)
            animal_class_name = sequence_dictionary[node.name][3]
            print(animal_class_name)
            if animal_class_name == 'Mammalia':
                node.set_style(mammalia)
            elif animal_class_name == 'Reptilia':
                node.set_style(reptilia)
            elif animal_class_name == 'Chondrichthyes':
                node.set_style(chondrichthyes)
            elif animal_class_name == 'Actinopterygii':
                node.set_style(actinopterygii)
            elif animal_class_name == 'Sarcopterygii':
                node.set_style(sarcopterygii)
            elif animal_class_name == 'Aves':
                node.set_style(aves)
            elif animal_class_name == 'Amphibia':
                node.set_style(amphibia)
            elif animal_class_name == 'Myxini':
                node.set_style(myxini)
            else:
                node.set_style(nstyle)
            # Set general leaf attributes
            node.set_style(general_leaf_style)
        else:
            node.set_style(nstyle)

    # ADD IMAGES
    #for key, value in dict_of_images.items():
    #    imgFace = ImgFace(IMG_PATH+value, height = 40)
    #    (t & key).add_face(imgFace, 0, "aligned")
    #    imgFace.margin_right = 10
    #    imgFace.hzalign = 2

    # ADD TEXT
    for key, value in sequence_dictionary.items():
        if key == outgroup_name:
            print(key, value[0], value[1], value[2], value[3])
            textFace = TextFace(value[1] + " (" + value[2] + ", " + key + ")   ", fsize = 16)
            (t & key).add_face(textFace, 2, "aligned")
            textFace.margin_left = 10
            print(key, value[0], value[1], value[2], value[3])
            textFace = TextFace(" ", fsize = 16)
            (t & key).add_face(textFace, 2, "aligned")
            textFace.margin_left = 10
            print(key, value[0], value[1], value[2], value[3])
            textFace = TextFace("Consensus score", fsize = 16)
            (t & key).add_face(textFace, 2, "aligned")
            textFace.margin_left = 10
        else:
            print(key, value[0], value[1], value[2], value[3])
            textFace = TextFace(value[1] + " (" + value[2] + ", " + key + ")   ", fsize = 16)
            (t & key).add_face(textFace, 2, "aligned")
            textFace.margin_left = 10

    # ADD DUMMY TEXT
    for key, value in sequence_dictionary.items():
        textFace3 = TextFace("   ", fsize = 16)
        (t & key).add_face(textFace3, 0, "aligned")

    # ADD SVG IMAGES
    for key, value in sequence_dictionary.items():
        svgFace = SVGFace(IMG_PATH+value[0], height = 40)
        (t & key).add_face(svgFace, 1, "aligned")
        svgFace.margin_right = 10
        svgFace.margin_left = 10
        svgFace.hzalign = 2
        animal_class_name = value[3]
        if animal_class_name == 'Mammalia':
            svgFace.background.color = "Chocolate"
        elif animal_class_name == 'Reptilia':
            svgFace.background.color = "Olive"
        elif animal_class_name == 'Chondrichthyes':
            svgFace.background.color = "SteelBlue"
        elif animal_class_name == 'Actinopterygii':
            svgFace.background.color = "CornflowerBlue"
        elif animal_class_name == 'Sarcopterygii':
            svgFace.background.color = "DarkCyan"
        elif animal_class_name == 'Aves':
            svgFace.background.color = "DarkSalmon"
        elif animal_class_name == 'Amphibia':
            svgFace.background.color = "DarkSeaGreen"
        elif animal_class_name == 'Myxini':
            svgFace.background.color = "LightBlue"
        else:
            svgFace.background.color = "White"


    color_dict = {  '-':"White",
                    '0':"#FF6666",
                    '1':"#EE7777",
                    '2':"#DD8888",
                    '3':"#CC9999",
                    '4':"#BBAAAA",
                    '5':"#AABBBB",
                    '6':"#99CCCC",
                    '7':"#88DDDD",
                    '8':"#77EEEE",
                    '9':"#66FFFF"}

    # ADD SEQUENCES AS TEXT (ADDING MORE THAN ONE SEQUENCE DOES NOT WORK)
    for key, value in dict_of_sequences.items():
        if key == outgroup_name:
            for char in range(0, len(value)):
                textFace2 = TextFace(value[char], fsize = 16)
                (t & key).add_face(textFace2, char+3, "aligned")
                textFace2.background.color = color_dict[dict_of_evaluations[key][char]]
                # This is for the consensus coloring, print a space instead of the sequence charactger of the outgroup sequence
                textFace4 = TextFace(' ', fsize = 16)
                (t & key).add_face(textFace4, char+3, "aligned")
                textFace4.background.color = "White"

                textFace5 = TextFace(' ', fsize = 16)
                (t & key).add_face(textFace5, char+3, "aligned")
                textFace5.background.color = color_dict[dict_of_evaluations['cons'][char]]

                textFace6 = TextFace(' ', fsize = 16)
                (t & key).add_face(textFace6, char+3, "aligned")
                textFace6.background.color = "White"
                # margins have the same color as the consensus color
                # These do not work!!!!
                #textFace4.margin_top = 10
                #textFace4.margin_bottom = 10

        else:
            for char in range(0, len(value)):
                #print(str(char) + " of " + str(len(value)))
                textFace2 = TextFace(value[char], fsize = 16)
                (t & key).add_face(textFace2, char+3, "aligned")
                #print('key: '+key)
                #print(dict_of_evaluations[key])
                #print(dict_of_evaluations[key][char])
                textFace2.background.color = color_dict[dict_of_evaluations[key][char]]

    # ROTATING SOME NODES CAN BE DONE HERE:
    # MOVE ACTINOPTERYGII NEXT TO THE OTHER FISH
    # Actinopterygii, e.g. Spotted gar - XP_006632034.2
    # Sarcopterygii, e.g. Coelacant - XP_006006690.1
    # Mammals, e.g. Rat - NP_446105.1
    # Birds, e.g. Chicken - XP_420532.3
    #n1 = t.get_common_ancestor("XP_420532.3", "NP_446105.1")

    # spotted gar and coelacant
    #n1 = t.get_common_ancestor("XP_006632034.2", "XP_006006690.1")
    #n1.swap_children()

    # Atlantic salmon and coelacant
    #n2 = t.get_common_ancestor("NP_001167218.1", "XP_006006690.1")
    #n2.swap_children()

    # Tibet frog and coelacant
    n2 = t.get_common_ancestor("XP_018419054.1", "XP_006006690.1")
    n2.swap_children()

    # Xenopus and coelacant
    n2 = t.get_common_ancestor("XP_002933363.1", "XP_006006690.1")
    n2.swap_children()

    # gekko and penguin
    n2 = t.get_common_ancestor("XP_015283812.1", "XP_009329004.1")
    n2.swap_children()

    # alligator and penguin
    n2 = t.get_common_ancestor("XP_006276984.1", "XP_009329004.1")
    n2.swap_children()


    # Add description to treefile
    #
    # Add fasta description line of outgroup sequence
    description_text = "Outgroup sequence: " + outgroup_name + "starfish VEGF-C\n"
    # Bootstrap analysis
    description_text += "Bootstrap: approximate Likelihood-Ratio Test\n"
    # Alignment methods
    description_text += "Alignment algorithm: m_coffee using" + "clustalw2, t_coffee, muscle, mafft, pcma, probcons"
    description_text += "\n"

    #ts.title.add_face(TextFace(description_text, fsize=12), column=0)
    t.render(SVG_TREEFILE, tree_style = ts, units = "mm", h = 120)

def run():
    if 'download' in args.options:
        download()
    if 'align' in args.options:
        align()
    if 'evaluate' in args.options:
        evaluate()
    if 'drawtree' in args.options:
        drawtree()

if __name__ == '__main__':
    run()
