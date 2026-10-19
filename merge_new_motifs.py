from pathlib import Path
import pandas as pd
from load_pdb_data import get_apfid_data


def segment_intersection_length(a1: int, a2: int, b1: int, b2: int) -> int:
    """
    Возвращает длину пересечения двух отрезков.

    a1, a2 - начало и конец первого отрезка
    b1, b2 - начало и конец второго отрезка
    """
    # Нормализуем отрезки, чтобы начало всегда было меньше конца
    a_start, a_end = min(a1, a2), max(a1, a2)
    b_start, b_end = min(b1, b2), max(b1, b2)

    # Начало пересечения — максимум из левых концов
    # Конец пересечения — минимум из правых концов
    intersect_start = max(a_start, b_start)
    intersect_end = min(a_end, b_end)

    # Если пересечение существует — возвращаем длину, иначе 0
    return max(0, intersect_end - intersect_start)

SOURCE_PATH = 'data/new_motifs'

NAMES = {
    '2.60.120 TNF-like': 'TNF-like',
    '1.20.1260 4-alpha bundle': '4α bundle',
    '3.30.370 Barnase': 'Barnase'
}

files = Path(SOURCE_PATH).rglob("*.csv")


allres = []

for file in files:
    print(file)
    df = pd.read_csv(file)
    df['motif_type'] = NAMES[file.parent.name]
    allres.append(df)

df = pd.concat(allres, ignore_index=True)
df.drop(columns=[
    'rmsd', 'fasta_score', 'stride_score', 'fasta_orig', 'stride_orig', 'stride_res', 'url'
], inplace=True)
df.sort_values(by=['motif_type', 'experiment_id', 'chain_id', 'start'], inplace=True)
df.drop_duplicates(['apfid', 'motif_type'])
df.rename(columns={
    'fasta_res': 'fasta',
    'start': 'motif_start',
    'end': 'motif_end',
}, inplace=True)
df['experiment_id'] = df['experiment_id'].apply(lambda x: x.upper())
df['source'] = 'PDB'
df['model_id'] = 0
df['full_chain'] = df['apfid'].apply(lambda x: len(x)==6)
df['motif_len'] = df['fasta'].apply(lambda x: len(x))
df.reset_index(inplace=True, drop=True)

results = []
for idx, row in df.iterrows():
    print(idx, row['motif_type'], row['apfid'])
    data = row.to_dict()
    if len(results) > 0:
        last_res = results[-1]
        if last_res['experiment_id'] == data['experiment_id'] and last_res['motif_type'] == data['motif_type'] and last_res['chain_id'] == data['chain_id']:
            intersect = segment_intersection_length(
                data['motif_start'],
                data['motif_end'],
                last_res['motif_start'],
                last_res['motif_end']
            )
            if intersect > 0 and intersect / min(last_res['motif_len'], data['motif_len']) > 0.3:
                continue
    apfid_data = get_apfid_data(data['apfid'])
    if apfid_data is None:
        continue
    data.update(apfid_data)
    results.append(data)

res_df = pd.json_normalize(results)
res_df.to_csv('data/DATA_NEW_MOTIFS.csv', index=False)


