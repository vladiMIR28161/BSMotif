import os
import numpy as np
import pandas as pd
import logomaker
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import io
import base64
import re
import seaborn as sns

def generate_heatmap_base64(score, branch_list):

    def extract_nums(s):
        nums = re.search(r'\{([^}]*)\}', s).group(1)
        parts = nums.split('.')
        parsed = []
        for p in parts:
            match = re.match(r'(\d+)([a-zA-Z]*)', p)
            number = int(match.group(1))
            suffix = match.group(2)
            parsed.append((number, suffix))
        return tuple(parsed)

    def extract_braces(s):
        match = re.search(r'\{[^}]*\}', s)
        return match.group(0) if match else s

    base_fontsize = 5
    pairs_list = []
    repeats = []
    branch_list = branch_list.split(", ")
    first_code_branch = re.findall(r"\{(.*?)\}", branch_list[0])

    if len(first_code_branch[0].split(".")) == 1:
        query = "Query_superclass"
        target = "Target_superclass"

    elif len(first_code_branch[0].split(".")) == 2:
        query = "Query_class"
        target = "Target_class"

    elif len(first_code_branch[0].split(".")) == 3:
        query = "Query_family"
        target = "Target_family"

    elif len(first_code_branch[0].split(".")) == 4:
        query = "Query_subfamily"
        target = "Target_subfamily"

    else:
        query = "Query_gene"
        target = "Target_gene"

    for branch1 in branch_list:
        for branch2 in branch_list:

            if [branch1, branch2] in repeats:
                continue

            pair_df = score.loc[
                (((score[query] == branch1) & (score[target] == branch2)) |
                 ((score[query] == branch2) & (score[target] == branch1)))
            ]

            med = np.median(pair_df.Score_TF.unique()).round(3)

            if np.isnan(med):
                med = '#'

            pairs_list.append([branch1, branch2, med])

            repeats.append([branch1, branch2])
            repeats.append([branch2, branch1])

    normalized = []

    for a, b, val in pairs_list:

        if extract_nums(a) > extract_nums(b):
            a, b = b, a
        normalized.append([a, b, val])
    pairs_list = sorted(normalized, key = lambda x: (extract_nums(x[0]), extract_nums(x[1])))
    df = pd.DataFrame(pairs_list, columns = ['row', 'col', 'val'])
    order = sorted(set(df['row']) | set(df['col']), key = extract_nums)
    df['row'] = pd.Categorical(df['row'], categories = order, ordered = True)
    df['col'] = pd.Categorical(df['col'], categories = order, ordered = True)
    result = df.pivot(index = 'row', columns = 'col', values = 'val')
    result = result.fillna('')
    result = result.rename_axis(index = None, columns = None)
    result.columns = result.columns.str.extract(r'(\{[^}]*\})', expand = False)
    heatmap_data = result.replace("#", np.nan).replace("", np.nan).astype(float)
    cell_size = 0.6
    n = len(order)
    fig, ax = plt.subplots(figsize = (n * cell_size, n * cell_size))
    cmap = plt.get_cmap('Reds', 100)

    sns.heatmap(heatmap_data, annot = False, vmin = 0, vmax = 3, cmap = cmap, linewidths = 0.8, linecolor = "grey",
        cbar_kws = {'shrink': 0.8, 'label': 'Median', 'pad': 0.05, 'location': "bottom"}, ax = ax)

    cbar = ax.collections[0].colorbar
    cbar.set_ticks([0, 1, 2, 3])
    cbar.set_ticklabels(["0", "1", "2", "3"])
    cbar.set_label('Median', fontsize = base_fontsize)
    cbar.ax.tick_params(labelsize = base_fontsize)

    for i in range(result.shape[0]):
        for j in range(result.shape[1]):
            value = result.iloc[i, j]
            if value != "":
                try:
                    num = float(value)
                    color = 'white' if num > 1.5 else 'black'
                except:
                    color = 'black'
                ax.text(j + 0.5, i + 0.5, str(value), ha = 'center', va = 'center', color = color, fontsize = base_fontsize)

    tick_fontsize = base_fontsize
    if result.shape[1] == 1:
        rotation_value = 0
    else:
        rotation_value = 90
    ax.set_xticklabels([extract_braces(x.get_text()) for x in ax.get_xticklabels()], rotation = rotation_value, ha = 'center', fontsize = tick_fontsize)
    ax.set_yticklabels([t.get_text() for t in ax.get_yticklabels()], rotation = 0, fontsize = tick_fontsize)
    ax.tick_params(axis = 'both', labelsize = tick_fontsize)
    ax.xaxis.tick_top()
    ax.set_aspect('equal')
    plt.subplots_adjust(left = 0.2, right = 0.9, top = 0.9, bottom = 0.2)
    buf = io.BytesIO()
    plt.savefig(buf, format = "png", dpi = 200, bbox_inches = 'tight', pad_inches = 0.05)
    buf.seek(0)
    b64 = base64.b64encode(buf.getvalue()).decode('utf-8')
    plt.close()
    return b64

def generate_motif_logos(results_df, classification_tsv, meme_dir, output_html, output_pcm_tmp, score_result, branch_result):

    def PCM_to_reverse_complement(PCM):
        reverse_complement_PCM = PCM[::-1][:, ::-1]
        return reverse_complement_PCM

    #Function for comparing all motifs in a dictionary with each other
    def compare_motifs(total_pcm, total_pcm2):
        results = {}

        for i in range(1, len(total_pcm) + 1):
            results[i] = {}

            for k in range(1, len(total_pcm2) + 1):
                results[i][k] = []

                for orientation in range(2):

                    fixed_motif = total_pcm[i][0]

                    if orientation == 0:
                        sliding_motif = total_pcm2[k][0]
                    else:
                        sliding_motif = PCM_to_reverse_complement(total_pcm2[k][0])

                    Lf = len(fixed_motif)
                    Ls = len(sliding_motif)

                    # the range of REAL shifts sliding can start from -(Ls-1) to Lf-1
                    for shift in range(-(Ls - 1), Lf):

                        score = 0
                        overlap = 0

                        # calculate the overlap
                        for j in range(Ls):
                            fi = j + shift
                            if 0 <= fi < Lf:
                                score += sum(abs(fixed_motif[fi] - sliding_motif[j]))
                                overlap += 1

                        # penalty for gaps
                        gap_penalty = (Ls - overlap) * 2

                        results[i][k].append(score + gap_penalty)

        return results

    # ALIGNMENT RESTORE FUNCTION
    def reconstruct_alignment(fixed, sliding, shift, bg=None):
        if bg is None:
            bg = np.array([0.25, 0.25, 0.25, 0.25])
        start = min(0, shift)
        end = max(len(fixed), shift + len(sliding))
        total_len = end - start
        aligned_fixed = []
        aligned_sliding = []
        for i in range(total_len):
            fi = i + start
            si = fi - shift
            if 0 <= fi < len(fixed):
                aligned_fixed.append(fixed[fi])
            else:
                aligned_fixed.append(bg)
            if 0 <= si < len(sliding):
                aligned_sliding.append(sliding[si])
            else:
                aligned_sliding.append(bg)
        return np.array(aligned_fixed), np.array(aligned_sliding)

    def merge_aligned_motifs(aligned_motif1, aligned_motif2):
        if aligned_motif1.shape != aligned_motif2.shape:
            raise ValueError("Мотивы должны быть одинаковой формы")
        merged_motif = (aligned_motif1 + aligned_motif2) / 2.0
        return merged_motif

    def normalize_pcm(pcm):
        pcm = pcm.copy()
        for i in range(len(pcm)):
            row_sum = pcm[i].sum()
            if row_sum > 0:
                pcm[i] /= row_sum
        return pcm

    # PCM to BITS (Information Content)
    def pcm_to_bits(pcm, background=None):
        if background is None:
            background = np.array([0.25, 0.25, 0.25, 0.25])
        bits = np.zeros_like(pcm)
        for i in range(len(pcm)):
            for j in range(4):
                if pcm[i, j] > 0:
                    bits[i, j] = pcm[i, j] * np.log2(pcm[i, j] / background[j])
                else:
                    bits[i, j] = 0.0
        return bits

    # Tail trimming
    def trim_by_information(bits_matrix, threshold=0.05):
        info_per_pos = bits_matrix.sum(axis=1)
        left = 0
        while left < len(info_per_pos) and info_per_pos[left] < threshold:
            left += 1
        right = len(info_per_pos) - 1
        while right >= 0 and info_per_pos[right] < threshold:
            right -= 1
        # If everything is cut off, return the original motif
        if left > right:
            return bits_matrix
        return bits_matrix[left:right + 1]

    # DATAFRAME FOR LOGOMAKER
    def bits_to_logomaker_df(bits_matrix):
        return pd.DataFrame(bits_matrix, columns=["A", "C", "G", "T"])

    branch = results_df.copy()
    branch["Logo_df"] = None
    branch["Heatmap_df"] = None

    file_path = classification_tsv
    df0 = pd.read_csv(file_path, sep='\t', index_col=0)

    level = ["Superclass", "Class", "Family", "Subfamily", "Gene"]

    for br in branch.List:
        print(br)
        smlrt = branch.loc[branch['List'] == br, 'Similarity'].values
        new_file = open (output_pcm_tmp , 'w')
        for i in range(len(br.split(', '))):
            for lvl in level:
                df1 = df0.loc[df0[lvl] == br.split(", ")[i]]
                if len(df1) != 0:
                    break
            for ind in df1.ID:
                flag = 0
                pwm = open(os.path.join(meme_dir, f"{ind}.meme"))
                for line in pwm:
                    if line.startswith("MOTIF "):
                        main_stroka = "> " + line.rstrip().split(" ")[1]
                    elif line.startswith("letter-probability"):
                        main_stroka += " length=" + line.split("w= ")[1].split(" nsites=")[0] + "\n"
                    elif line.startswith("0.") or line.startswith("1.") or line.startswith(" 0.") or line.startswith(" 1."):
                        if flag == 0:
                            new_file.write(main_stroka)
                            new_file.write(line.lstrip())
                            flag = 1
                        else:
                            new_file.write(line.lstrip())
                new_file.write("\n")
        new_file.close()


        new_file = open (output_pcm_tmp)
        motif_and_matrix = {}
        N = 0
        for line in new_file:
            if line.startswith("> "):
                N += 1
                key = line.split()[1]
                motif_and_matrix[N] = []
            elif line.startswith("0.") or line.startswith("1."):
                values = list(map(float, line.split()))
                motif_and_matrix[N].append(values)
        new_file.close()

        for key in motif_and_matrix:
            motif_and_matrix[key] = np.array(motif_and_matrix[key])

        background = [0.25, 0.25, 0.25, 0.25]

        first_key, first_value = next(iter(motif_and_matrix.items()))
        lst1 = []
        lst1.append(first_value)

        for key, value in motif_and_matrix.items():
            motif_1 = {1: lst1}
            lst = []
            lst.append(value)
            motif_2 = {1: lst}
            if len(motif_2[1][0]) > len(motif_1[1][0]):
                motif1 = motif_2
                motif2 = motif_1
            else:
                motif1 = motif_1
                motif2 = motif_2
            results = compare_motifs(motif1, motif2)
            print(motif1)
            print(motif2)
            print(results)

            # FINDING THE BEST ALIGNMENT
            i = 1
            k = 1
            scores = results[i][k]
            best_index = int(np.argmin(scores))
            best_score = scores[best_index]
            fixed = motif1[i][0]
            sliding = motif2[k][0]
            num_shifts = len(fixed) + len(sliding) - 1
            orientation = best_index // num_shifts
            shift_index = best_index % num_shifts
            best_shift = shift_index - (len(sliding) - 1)
            fixed = motif1[i][0]

            if orientation == 0:
                sliding = motif2[k][0]
            else:
                sliding = PCM_to_reverse_complement(motif2[k][0])

            aligned_motif1, aligned_motif2 = reconstruct_alignment(fixed, sliding, best_shift)

            print("Best score:", best_score)
            print("Best shift:", best_shift)
            print("Aligned motif 1 shape:", aligned_motif1.shape)
            print("Aligned motif 2 shape:", aligned_motif2.shape)

            merged_motif = merge_aligned_motifs(aligned_motif1, aligned_motif2)
            print("Merged motif shape:", merged_motif.shape)
            print(merged_motif)
            merged_motif_normalized = normalize_pcm(merged_motif)
            bits_matrix = pcm_to_bits(merged_motif_normalized)
            bits_matrix[bits_matrix < 0] = 0
            trimmed_bits = trim_by_information(bits_matrix, threshold=0.25)
            logo_df = bits_to_logomaker_df(trimmed_bits)

            logo_df.index = np.arange(1, len(logo_df) + 1)
            lst1.clear()
            lst1 = []
            lst1.append(merged_motif)

        # LOGO
        plt.figure(figsize=(len(logo_df) * 0.6, 4))
        logomaker.Logo(logo_df)

        ax = plt.gca()
        ax.set_xticks(logo_df.index)
        ax.set_xticklabels(logo_df.index)

        ax.set_ylim(0, 2)

        plt.ylabel("Information (bits)")
        plt.xlabel("Position")
        plt.title("Merged motif")
        plt.tight_layout()
        buf = io.BytesIO()
        plt.savefig(buf, format="png", dpi=300)
        buf.seek(0)

        b64_string = base64.b64encode(buf.getvalue()).decode('utf-8')
        branch.loc[branch['List'] == br, 'Logo_df'] = b64_string
        heatmap_b64 = generate_heatmap_base64(score_result, br)
        branch.loc[branch['List'] == br, 'Heatmap_df'] = heatmap_b64
        plt.close()

    os.remove(output_pcm_tmp)
    html_content = """
    <html>
    <head>
    <meta charset="utf-8">
    <title>Table of branches</title>

    <style>

    body {
        font-family: Arial, sans-serif;
        margin: 20px;
    }

    table {
        border-collapse: collapse;
        width: 100%;
        table-layout: fixed;
    }

    th, td {
        border: 1px solid black;
        padding: 10px;
        vertical-align: middle;
        text-align: center;
    }

    th {
        background-color: #f2f2f2;
        font-size: 20px;
    }

    td {
        font-size: 18px;
    }

    th:nth-child(1), td:nth-child(1) {
        width: 7%;
        white-space: nowrap;
    }

    th:nth-child(2), td:nth-child(2) {
        width: 15%;
        text-align: center;
    }

    th:nth-child(3), td:nth-child(3) {
        width: 7%;
        white-space: nowrap;
        font-size: 20px;
    }

    th:nth-child(4), td:nth-child(4) {
        width: 36%;
    }

    th:nth-child(5), td:nth-child(5) {
        width: 35%;
    }

    img {
        display: block;
        margin-left: auto;
        margin-right: auto;
        max-width: 100%;
        height: auto;
    }

    </style>

    </head>
    <body>
    <h1>Table of branches of motifs</h1>
    <table>
    <tr>
    <th>Branch</th>
    <th>List</th>
    <th>Similarity</th>
    <th>Logo</th>
    <th>Heatmap</th>
    </tr>
    """

    for _, row in branch.iterrows():

        if pd.isna(row['Logo_df']):
            continue

        html_content += f"""
        <tr>
            <td>{row['Branch']}</td>
            <td>{row['List']}</td>
            <td>{row['Similarity']}</td>
            <td>
                <img src="data:image/png;base64,{row['Logo_df']}">
            </td>
            <td>
                <img src="data:image/png;base64,{row['Heatmap_df']}">
            </td>
        </tr>
        """

    html_content += """
    </table>
    </body>
    </html>
    """

    with open(output_html, "w", encoding="utf-8") as f:
        f.write(html_content)
