
""" This script calculates the features needed for the classifier. It uses the gene coordinate file from prodigal
 or a .gff file. """

import pandas as pd
import numpy as np
import os


def add_features(dictionary, genelist, filetype, seqlength):
    data_dict = dictionary
    gene_list = genelist

    genes_same = 0  # <- <- or -> -> (++, --)
    genes_diff = 0  # <- -> (-+)
    genes_general = []  # All orientations together
    data_dict["nr_genes"].append(len(gene_list))

    if len(gene_list) < 2:
        data_dict["ratio_same_orientation"].append(np.nan)

        data_dict["ID_general_avg"].append(np.nan)
        data_dict["ID_general_std"].append(np.nan)
        data_dict["ID_Q1"].append(np.nan)
        data_dict["ID_Q3"].append(np.nan)

    if len(gene_list) >= 2:

        for i in range(0, (len(gene_list) - 1)):
            gene = gene_list[i]
            next_gene = gene_list[i + 1]
            distance = ''
            gene_orientation = ''
            next_gene_orientation = ''

            if filetype == "coord":
                distance = int(next_gene[1]) - int(gene[2])
                gene_orientation = gene[3]
                next_gene_orientation = next_gene[3]

            if filetype == "GFF":
                distance = int(next_gene[0]) - int(gene[1])
                gene_orientation = gene[2]
                next_gene_orientation = next_gene[2]

            if distance < 0:
                distance = 0

            genes_general.append(distance)
            orientation = str(gene_orientation) + "_" + str(next_gene_orientation)

            if orientation == "+_+" or orientation == "-_-":
                genes_same += 1

            if orientation == "-_+" or orientation == "+_-":
                genes_diff += 1

        ratio_same = genes_same / (genes_same + genes_diff)
        data_dict["ratio_same_orientation"].append(ratio_same)

        if len(genes_general) > 1:
            data_dict["ID_general_avg"].append(np.mean(genes_general))
            data_dict["ID_general_std"].append(np.std(genes_general))
            data_dict["ID_Q1"].append(np.quantile(genes_general, 0.25))
            data_dict["ID_Q3"].append(np.quantile(genes_general, 0.75))

        if len(genes_general) == 1:
            data_dict["ID_general_avg"].append(genes_general[0])
            data_dict["ID_general_std"].append(0)
            data_dict["ID_Q1"].append(0)
            data_dict["ID_Q3"].append(0)

        if len(genes_general) < 1:
            data_dict["ID_general_avg"].append(np.nan)
            data_dict["ID_general_std"].append(np.nan)
            data_dict["ID_Q1"].append(np.nan)
            data_dict["ID_Q3"].append(np.nan)

    length = []

    if filetype == "GFF":
        for gene in gene_list:
            gene_length = int(gene[1]) - int(gene[0])
            length.append(gene_length)

    if filetype == "coord":
        for gene in gene_list:
            gene_length = int(gene[2]) - int(gene[1])
            length.append(gene_length)

    density = sum(length) / seqlength
    data_dict["gene_density"].append(density)

    if len(length) != 0:
        avglength = sum(length) / len(length)
        data_dict["gene_length"].append(avglength)

    if len(length) == 0:
        avglength = sum(length)
        data_dict["gene_length"].append(avglength)


def calc_features(contig_file, outfile):

    data_dict = {"contig": [],
                 "contig_length": [],
                 "nr_genes": [],
                 "ratio_same_orientation": [],
                 "ID_general_avg": [],
                 "ID_general_std": [],
                 "gene_density": [],
                 "gene_length": [],
                 "ID_Q1": [],
                 "ID_Q3": [],
                 "rbs_ratio": []}

    gene_list = "empty"
    rbs_list = "empty"
    filetype = "unknown"

    if ".gff" in contig_file:
        filetype = "GFF"
        print("File extension: GFF")

    if ".gene" in contig_file:
        filetype = "coord"
        print("File extension: .genes")

    with open(contig_file, newline='') as coords_file:
        for line in coords_file:

            if line.startswith("# Seq"):

                if gene_list != "empty":   # This is the gene list from the previous contig
                    if len(rbs_list) == 0:
                        rbs_ratio = np.nan
                        print("RBS_list = 0", rbs_ratio)
                    if len(rbs_list) > 0:
                        rbs_ratio = 1 - (rbs_list.count("None") / len(rbs_list))

                    data_dict["rbs_ratio"].append(rbs_ratio)

                    add_features(data_dict, gene_list, filetype, seqlength)

                header = line
                seqname = header.split(";")[2].split('"')[1].split(" ")[0]
                seqlength = int(header.split(";")[1].split("=")[1])

                data_dict["contig"].append(seqname)
                data_dict["contig_length"].append(seqlength)

                gene_list = []

            if filetype == "coord":
                if line.startswith(">"):
                    info = line.split("\n")[0].split(">")[1].split("_")
                    gene_list.append(info)

            if filetype == "GFF":
                if not line.startswith("#"):
                    stats = line.split("\t")[8]
                    rbs_motifs = stats.split(";")[3].split("=")[1]
                    rbs_list.append(rbs_motifs)

                    info = [line.split("\t")[3], line.split("\t")[4], line.split("\t")[6]]
                    gene_list.append(info)

        # This is for the last contig:
        if gene_list != "empty":
            if len(rbs_list) == 0:
                rbs_ratio = np.nan
            if len(rbs_list) > 0:
                rbs_ratio = 1 - (rbs_list.count("None") / len(rbs_list))

            data_dict["rbs_ratio"].append(rbs_ratio)
            add_features(data_dict, gene_list, filetype, seqlength)

    # Saving the dictionary to a dataframe that will be used for the classifier:
    df = pd.DataFrame(data_dict, columns=list(data_dict.keys()))

    df.to_csv(os.path.join(outfile, "featuretable.csv"), index=False)
