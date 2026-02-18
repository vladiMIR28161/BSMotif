import os
import numpy as np
import pandas as pd
import logomaker
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import io
import base64

def generate_motif_logos(results_df, classification_tsv, meme_dir, output_html, output_pcm_tmp):

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
        plt.close()

    os.remove(output_pcm_tmp)
    html_content = """
    <html>
    <head>
    <meta charset="utf-8">
    <title>Table of branches</title>

    <style>
    table {
        border-collapse: collapse;
        width: 100%;
        table-layout: fixed;   /* фиксированная ширина колонок */
    }

    th, td {
        border: 1px solid black;
        padding: 5px;
        vertical-align: middle;
    }

    th:nth-child(1), td:nth-child(1) {
        width: 6%;
        white-space: nowrap;   /* запрет переноса */
    }

    th:nth-child(2), td:nth-child(2) {
        width: 59%;
    }

    th:nth-child(3), td:nth-child(3) {
        width: 5%;
        text-align: center;
    }

    th:nth-child(4), td:nth-child(4) {
        width: 30%;
        text-align: center;
    }

    img {
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
        </tr>
        """

    html_content += """
    </table>
    </body>
    </html>
    """

    with open(output_html, "w", encoding="utf-8") as f:
        f.write(html_content)
