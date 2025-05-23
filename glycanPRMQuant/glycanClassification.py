import pandas as pd

def classifyGlycan(consolidated_csv: str) -> pd.DataFrame:
    """
    Read a consolidated AUC CSV with a 'Glycan' column of 5-digit integer
    compositions, split that into five monosaccharide counts, and assign
    each glycan to one of:
      - high mannose
      - sialofucosylated
      - fucosylated
      - sialylated
      - other

    Returns the original DataFrame plus columns pos1…pos5 and 'Class'.
    """
    df = pd.read_csv(consolidated_csv)

    # ensure Glycan is a 5-digit string, pad with leading zeros if needed
    df['gly_str'] = df['Glycan'].astype(str).str.zfill(5)

    # split into five integer columns: pos1…pos5
    for i in range(5):
        df[f'pos{i+1}'] = df['gly_str'].str[i].astype(int)

    # classification function
    def _classify(row):
        # high mannose: HexNAc2 (pos1==2) & Hexose>=5 (pos2>=5)
        if row['pos1'] == 2 and row['pos2'] >= 5:
            return 'high mannose'
        # sialofucosylated: both fucose (pos4) and sialic acid (pos5)
        if row['pos4'] > 0 and row['pos5'] > 0:
            return 'sialofucosylated'
        # fucosylated only
        if row['pos4'] > 0:
            return 'fucosylated'
        # sialylated only
        if row['pos5'] > 0:
            return 'sialylated'
        return 'other'

    df['Class'] = df.apply(_classify, axis=1)

    # drop helper string column if you like
    df = df.drop(columns=['gly_str'])

    return df
