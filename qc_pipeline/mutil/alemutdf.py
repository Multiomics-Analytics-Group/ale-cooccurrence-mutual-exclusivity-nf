import os
import os.path
import re
import pandas as pd
import numpy as np

from mutil.mut import get_del_size, get_ins_size, get_inv_size, get_con_size, get_sub_size, get_amp_size, is_frameshift, is_coding_mut, is_premature_stop_codon_SNP, is_readthrough_codon_SNP, is_start_codon_removal


# TODO: Integrate functionality to get the series of fixed mutations defined by a final mutation with a frequency larger than a given threshold into this list of scripts. This functionality is implemented in substitution_rate_decline.ipynb.


def get_coding_muts_df(mut_df):
    coding_mut_df = mut_df[~mut_df["Details"].str.contains("intergenic")
                           & ~mut_df["Details"].str.contains("pseudogene")
                           & ~mut_df["Details"].str.contains("noncoding")].copy()
    return coding_mut_df


def get_noncoding_muts_df(mut_df):
    coding_mut_df = mut_df[mut_df["Details"].str.contains("intergenic")
                           | mut_df["Details"].str.contains("pseudogene")
                           | mut_df["Details"].str.contains("noncoding")].copy()
    return coding_mut_df


def get_genetic_muts_df(mut_df):
    df = mut_df[~mut_df["Details"].str.contains("intergenic")].copy()
    return df


def get_mut_dataframe(CSV_file_path,
                      log,
                      include_dups = False,
                      intragenic_muts_only = False):
    """Mutations of one ALEdb file. Returns (mutations, a list of the removed rows with their
    removal_reason, (number of sample columns dropped for repeating another column's name, how many
    of them differ from that column, how many calls those hold)). log(message) records each note
    about the file."""

    # Step 1: Import database
    try:
        raw_db = pd.read_csv(CSV_file_path)
    except pd.errors.EmptyDataError:
        log(f"{CSV_file_path} is empty, skipped")
        return pd.DataFrame(), [], (0, 0, 0)
    
    # Check if the dataframe is empty after reading
    if raw_db.empty or len(raw_db) == 0:
        log(f"{CSV_file_path} contains no data, skipped")
        return pd.DataFrame(), [], (0, 0, 0)

    # Clean duplicate columns (remove .1, .2, .3, etc. suffixes). pandas adds the suffix to a column
    # whose name repeats an earlier one. Each dropped column is compared with the column it repeats,
    # so the log says when it held other calls.
    is_dup = raw_db.columns.str.match(r'.*\.\d+$')
    dup_cols = list(raw_db.columns[is_dup])
    differing = [c for c in dup_cols if not raw_db[c].equals(raw_db.get(re.sub(r'\.\d+$', '', c)))]
    lost_calls = int(sum(raw_db[c].notnull().sum() for c in differing))
    raw_db = raw_db.loc[:, ~is_dup]
    if dup_cols:
        log(f"removed {len(dup_cols)} sample columns whose names repeat another column (name.1) "
            f"from {CSV_file_path}: {len(dup_cols) - len(differing)} identical to the column they "
            f"repeat, {len(differing)} different, with {lost_calls} calls that are not loaded")

    if 'Function' in raw_db.columns:
        raw_db = raw_db.drop('Function', axis=1)
    if 'Product' in raw_db.columns:
        raw_db = raw_db.drop('Product', axis=1)
    if 'GO Process' in raw_db.columns:
        raw_db = raw_db.drop('GO Process', axis=1)
    if 'GO Component' in raw_db.columns:
        raw_db = raw_db.drop('GO Component', axis=1)
    
    # TODO: The below will crash if most than these columns. Needs to ignore other mutation columns.
    # Step 2: Separate columns based on usage
    keep_cols = ['Position','Mutation Type','Sequence Change','Details','Gene (Scrollable)']
    if  "Reference Seq" in raw_db.columns:  # For backwards compatibility with older exported mutation data
        keep_cols.append("Reference Seq")
    if "Mut ID" in raw_db.columns:  # For backwards compatibility with older exported mutation data
        keep_cols.append("Mut ID")
    
    # Check if required columns exist
    missing_cols = [col for col in ['Position','Mutation Type','Sequence Change','Details','Gene (Scrollable)'] if col not in raw_db.columns]
    if missing_cols:
        raise ValueError(f"{CSV_file_path} is missing required columns: {missing_cols}")

    mut_cols = sorted(list(set(raw_db.columns) - set(keep_cols)))

    # Check if there are any mutation columns
    if not mut_cols:
        raise ValueError(f"{CSV_file_path} contains no mutation columns")
    
    file_name = os.path.basename(CSV_file_path)
    exp_name_from_file = file_name.replace(".csv", '')
    
    # Step 3: Shift mutation column names into row identifiers
    csv_file_mutat_df = pd.DataFrame()
    for col in mut_cols:
        df = raw_db[raw_db[col].notnull()][keep_cols]
        
        # Skip if no non-null values in this column
        if df.empty:
            continue
            
        exp_name = '_'.join(col.split(' ')[:-4])
        if exp_name == '':  # Will happen with newer mutation data exported from ALEdb
            exp_name = exp_name_from_file
        df['exp'] = exp_name
        
        # Handle potential issues with column name parsing
        try:
            df['ale'] = int(col.split(' ')[-4][1:])
            df['flask'] = int(col.split(' ')[-3][1:])
            df['isolate'] = int(col.split(' ')[-2][1:])
            df['tech_rep'] = int(col.split(' ')[-1][1:])
        except (IndexError, ValueError) as e:
            raise ValueError(f"could not parse column name '{col}' in {CSV_file_path}") from e
            
        df['presence'] = raw_db[raw_db[col].notnull()][col]
        csv_file_mutat_df = pd.concat([csv_file_mutat_df,df])

    # A file with mutation rows but no call in any sample column is not an empty deposit
    if csv_file_mutat_df.empty:
        raise ValueError(f"{CSV_file_path} has mutation rows but no value in any sample column")

    csv_file_mutat_df = csv_file_mutat_df[['exp','ale','flask','isolate','tech_rep','presence'] + keep_cols]
    csv_file_mutat_df = csv_file_mutat_df.fillna('')
    
    # Remove mutation entries with empty gene since they will screw up mutat_df.groupby(['Gene', ...])
    no_gene = csv_file_mutat_df['Gene (Scrollable)'] == ''
    removed = [csv_file_mutat_df.loc[no_gene].assign(removal_reason="no gene name")]
    csv_file_mutat_df = csv_file_mutat_df.loc[~no_gene]
    
    # Remove weird characters between gene names in multiple gene annotation.
    csv_file_mutat_df['Gene (Scrollable)'] = csv_file_mutat_df['Gene (Scrollable)'].str.replace("  ", " ")
    
    if not include_dups:
        dup = csv_file_mutat_df['Details'] == 'Duplication'
        removed.append(csv_file_mutat_df.loc[dup].assign(removal_reason="Details is Duplication"))
        csv_file_mutat_df = csv_file_mutat_df.loc[~dup]
    
    if intragenic_muts_only:
        csv_file_mutat_df = csv_file_mutat_df.loc[csv_file_mutat_df['Gene (Scrollable)'].str.contains(',') == False]
    
    return csv_file_mutat_df, [r for r in removed if len(r)], (len(dup_cols), len(differing), lost_calls)


def get_all_sample_mut_df(dir_path,
                          log,
                          file_names,
                          include_dups = False,
                          intragenic_muts_only = False):
    """Mutations of the ALEdb files in dir_path whose names are in file_names. Returns (mutations, the
    rows removed while loading with removal_reason, counts of the files). Every row, loaded or removed,
    holds the name of its file in mut_file. log(message) records each note, such as a loaded or
    skipped file or dropped duplicate columns."""
    mutat_df = pd.DataFrame()
    mutat_df_list = []
    removed_list = []
    counts = {"files_loaded": 0, "files_without_rows": 0, "files_with_duplicated_columns": 0,
              "duplicated_columns_removed": 0, "duplicated_columns_differing": 0,
              "calls_in_differing_duplicated_columns": 0}

    # Get list of CSV files in sorted order of their names (by character code), so the row
    # order does not depend on the order in which the file system lists them
    csv_files = sorted(f for f in os.listdir(dir_path) if f.endswith('.csv'))
    # Load only the given file names
    csv_files = [f for f in csv_files if f in file_names]

    if not csv_files:
        raise ValueError(f"No CSV files found in {dir_path} among the selected file names")
    
    for file_name in csv_files:
        file_path = os.path.join(dir_path, file_name)
        log(f"loading {file_path}")
        
        # Get mutation dataframe, which will return empty DataFrame if file is empty
        mut_df, removed, (n_dup, n_differing, n_lost) = get_mut_dataframe(
            file_path, log, include_dups, intragenic_muts_only)
        counts["files_with_duplicated_columns"] += bool(n_dup)
        counts["duplicated_columns_removed"] += n_dup
        counts["duplicated_columns_differing"] += n_differing
        counts["calls_in_differing_duplicated_columns"] += n_lost
        # Two files can hold the same experiment name, so each row keeps the name of its file.
        removed_list.extend(r.assign(mut_file=file_name) for r in removed)
        
        # Only append if the dataframe is not empty
        if not mut_df.empty:
            mutat_df_list.append(mut_df.assign(mut_file=file_name))
            counts["files_loaded"] += 1
        else:
            log(f"skipped {file_name}: no mutation rows to load")
            counts["files_without_rows"] += 1
    
    # Concatenate all non-empty dataframes
    if not mutat_df_list:
        raise ValueError(f"No mutation data found in any CSV file in {dir_path}")
    mutat_df = pd.concat(mutat_df_list, ignore_index=True)

    if removed_list:
        removed_df = pd.concat(removed_list, ignore_index=True)
    else:
        # No row was removed. The columns are those of a removed row, with removal_reason before mut_file.
        removed_df = pd.DataFrame(columns=list(mutat_df.columns.drop("mut_file")) + ["removal_reason", "mut_file"])
    return mutat_df, removed_df, counts


def _get_exp_ale_set(mut_df):
    exp_ale_df = mut_df.copy()
    exp_ale_df["exp ale"] = exp_ale_df["exp"] + ' ' + exp_ale_df["ale"].map(str)
    exp_ale_set = set(exp_ale_df['exp ale'].tolist())
    return exp_ale_set


def get_gene_mut_mat(exp_ale_mut_gene_df):
    mut_df = exp_ale_mut_gene_df.copy()
    column_to_delete_list = ["presence",
                         "tech_rep",
                         "isolate",
                         "flask",
                         "Position",
                         "Mutation Type",
                         "Sequence Change",
                         "Details"]
    current_columns = list(mut_df.columns.values)
    for column_to_delete in column_to_delete_list:
        if column_to_delete in current_columns:
            del mut_df[column_to_delete]
    mut_df = mut_df.drop_duplicates()
    mut_df["exp ale"] = mut_df["exp"]+' '+mut_df["ale"].map(str)
    
    mut_df_column_name_set = _get_exp_ale_set(mut_df)

    # Get mut_mat_df indexes of Gene names.
    # Can't simply use the "Gene (Scrollable)" column of trunc_mut_df since a gene may be mutated in
    # more than one ALE exp, and we want a set of unique gene names.
    mut_df_index_set = set(mut_df["Gene (Scrollable)"].tolist())

    mut_mat_df = pd.DataFrame(columns=mut_df_column_name_set, index=mut_df_index_set)
    mut_mat_df = mut_mat_df.fillna(0)

    for gene_name, all_gene_mut_df in mut_df.groupby("Gene (Scrollable)"):
        for index, mut_df_row in all_gene_mut_df.iterrows():
            mut_mat_df.loc[gene_name, mut_df_row["exp ale"]] = 1
    
    return mut_mat_df


def get_gene_mut_count_mat(exp_ale_mut_gene_df):
    mut_df = exp_ale_mut_gene_df.copy()
    mut_df["exp ale"] = mut_df["exp"]+' '+mut_df["ale"].map(str)
    
    mut_df_column_name_set = _get_exp_ale_set(mut_df)

    # Get mut_mat_df indexes of Gene names.
    # Can't simply use the "Gene (Scrollable)" column of trunc_mut_df since a gene may be mutated in
    # more than one ALE exp, and we want a set of unique gene names.
    mut_df_index_set = set(mut_df["Gene (Scrollable)"].tolist())

    mut_mat_df = pd.DataFrame(columns=mut_df_column_name_set, index=mut_df_index_set)
    mut_mat_df = mut_mat_df.fillna(0)

    for gene_name, all_gene_mut_df in mut_df.groupby("Gene (Scrollable)"):
        for index, mut_df_row in all_gene_mut_df.iterrows():
            mut_mat_df.loc[gene_name, mut_df_row["exp ale"]] += 1
    
    return mut_mat_df


def get_mut_mat(exp_ale_mut_gene_df):
    mut_df = exp_ale_mut_gene_df.copy()
    column_to_delete_list = ["presence",
                         "tech_rep",
                         "isolate",
                         "flask",
                         "Position",
                         "Mutation Type",
                         "Details"]
    current_columns = list(mut_df.columns.values)
    for column_to_delete in column_to_delete_list:
        if column_to_delete in current_columns:
            del mut_df[column_to_delete]
    mut_df = mut_df.drop_duplicates()
    mut_df["exp ale"] = mut_df["exp"] + ' ' + mut_df["ale"].map(str)
    mut_df["gene seq change"] = mut_df["Gene (Scrollable)"] + ' ' + mut_df["Sequence Change"]
    
    mut_df_column_name_set = _get_exp_ale_set(mut_df)  # Get unique set of exp+ALE#

    mut_df_index_set = set(mut_df["gene seq change"].tolist())  # Get unique set of gene+(seq change) names.
    mut_mat_df = pd.DataFrame(columns=mut_df_column_name_set, index=mut_df_index_set)
    mut_mat_df = mut_mat_df.fillna(0)

    for mut, all_mut_df in mut_df.groupby("gene seq change"):
        for index, mut_df_row in all_mut_df.iterrows():
            mut_mat_df.loc[mut, mut_df_row["exp ale"]] = 1
    
    return mut_mat_df


def get_enrichment_muts(mut_df):
    trunc_mut_df = mut_df.copy()

    # If we are going to keep Duplications, though we want to remove the '[' and ']' from their gene annotations.
#     trunc_mut_df["Gene (Scrollable)"] = trunc_mut_df["Gene (Scrollable)"].map(lambda x: x.lstrip('[').rstrip(']'))

    # Removing duplications
    trunc_mut_df = trunc_mut_df[trunc_mut_df["Details"] != "Duplication"]

    # Removing unused columns
    del trunc_mut_df["tech_rep"]
    del trunc_mut_df["isolate"]
    # Could have the same mutation, but with a different presence due to differences between clonal and population reseq'ing
    del trunc_mut_df["presence"]
    trunc_mut_df = trunc_mut_df.drop_duplicates()

    dfs = []  # Create a list to store DataFrames
    for gene_mut_groupby in trunc_mut_df.groupby(["exp", "Gene (Scrollable)"]):
        gene_mut_df = gene_mut_groupby[1]
        mutation_count = gene_mut_df.shape[0]
        if mutation_count > 1:
            dfs.append(gene_mut_df)

    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()


def get_ALE_final_flask_df(ALE_df):
    final_flask_num = ALE_df["flask"].max()
    final_flask_df = ALE_df[ALE_df["flask"]==final_flask_num]
    return final_flask_df


def get_exp_final_flask_df(exp_df):
    dfs = []  # Create a list to store DataFrames
    for ale, ale_mut_df in exp_df.groupby("ale"):
        ALE_final_flask_df = get_ALE_final_flask_df(ale_mut_df)
        dfs.append(ALE_final_flask_df)
    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()


def _get_max_freq_mut_df(mut_df):
    max_mut_freq = mut_df["presence"].max()
    max_mut_df = mut_df[mut_df["presence"] == max_mut_freq]
    max_mut_df = max_mut_df.sort_values("flask")
    earliest_flask_max_mut_df = max_mut_df[:1]
    return earliest_flask_max_mut_df  # only need to return first row of max_mut_df


def get_ALE_max_freq_mut_df(ALE_mut_df, endpoint_flask_only, representative_isolate_only=False):
    if endpoint_flask_only:
        ALE_mut_df = get_ALE_final_flask_df(ALE_mut_df)

    # Select representative isolate (most unique mutations)
    if representative_isolate_only and len(ALE_mut_df) > 0:
        mutation_cols = ["Position", "Mutation Type", "Sequence Change", "Details", "Gene (Scrollable)"]
        
        unique_counts = ALE_mut_df.groupby('isolate').apply(
            lambda x: len(x[mutation_cols].drop_duplicates()),
            include_groups=False
        )
        
        best_isolate = unique_counts.idxmax()
        ALE_mut_df = ALE_mut_df[ALE_mut_df['isolate'] == best_isolate]
    
    dfs = []  # Create a list to store DataFrames
    #  List of everything that defines a mutation besides instance values (freq, isolate#, etc.)
    mutation_descriptors_list = ["Position", "Mutation Type", "Sequence Change", "Details", "Gene (Scrollable)"]
    
    for mut_group in ALE_mut_df.groupby(mutation_descriptors_list): 
        mut_df = mut_group[1]
        if len(mut_df) > 0:
            dfs.append(_get_max_freq_mut_df(mut_df))
    
    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()


def get_exp_max_freq_mut_df(exp_mut_df, endpoint_flask_only, representative_isolate_only=False):
    dfs = []
    for ale_name, ALE_mut_df in exp_mut_df.groupby(["ale"]):
        if len(ALE_mut_df) > 0:
            dfs.append(get_ALE_max_freq_mut_df(ALE_mut_df, endpoint_flask_only, representative_isolate_only))
    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()


def get_multi_exp_max_freq_mut_df(mut_df, endpoint_flask_only, representative_isolate_only=False):
    dfs = []  # Create a list to store DataFrames
    for group_name, exp_ale_mut_df in mut_df.groupby(["exp", "ale"]):
        if len(exp_ale_mut_df) > 0:
            dfs.append(get_ALE_max_freq_mut_df(exp_ale_mut_df, endpoint_flask_only, representative_isolate_only))
    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()


def get_filtered_mut_df(fixed_mut_df, filter_mut_df):
    if "flask" in filter_mut_df.columns:
        filter_mut_df = filter_mut_df.drop("flask", axis=1)
    if "isolate" in filter_mut_df.columns:
        filter_mut_df = filter_mut_df.drop("isolate", axis=1)
    if "tech_rep" in filter_mut_df.columns:
        filter_mut_df = filter_mut_df.drop("tech_rep", axis=1)
    if "presence" in filter_mut_df.columns:
        filter_mut_df = filter_mut_df.drop("presence", axis=1)
    return pd.merge(fixed_mut_df,
                    filter_mut_df,
                    how="inner",
                    on=["exp", "ale", "Position", "Mutation Type", "Sequence Change", "Details", "Gene (Scrollable)"])


# Filtering for fixed mutation SERIES that have a max freq larger than mut_freq_floor.
def get_filtered_fixed_mut_series_df(fixed_mut_df, mut_freq_floor):
    
    # TODO: Seems like building the filter_mut_df could be replaced with _get_max_freq_mut_df(...),
    # though need to make unit test to clarify.
    filter_mut_df_list = []
    for exp_name, exp_mut_df in fixed_mut_df.groupby(["exp"]):
        filter_mut_df_list.append(get_exp_max_freq_mut_df(exp_mut_df, endpoint_flask_only=True, representative_isolate_only=False))
    filter_mut_df = pd.concat(filter_mut_df_list)
    
    filter_mut_df = filter_mut_df[filter_mut_df["presence"]>mut_freq_floor]
    return get_filtered_mut_df(fixed_mut_df, filter_mut_df)

'''
# !!! get_ALE_mut_type_frac_df AND get_exp_avg_ale_mut_type_frac_df ARE OBSOLETE.
# USE get_mut_type_avg_frac_across_class_df INSTEAD

def get_ALE_mut_type_frac_df(mut_df, df_column_name):
    mut_type_set = set(mut_df[df_column_name])
    mut_type_frac_df = pd.DataFrame()
    for exp_ale_tuple, exp_ale_max_mut_df in mut_df.groupby(["exp", "ale"]):
        for mut_type in mut_type_set:
            mut_type_count = len(exp_ale_max_mut_df[exp_ale_max_mut_df[df_column_name]==mut_type])
            mut_type_fraction = mut_type_count/len(exp_ale_max_mut_df)

            mut_type_frac_df = mut_type_frac_df.append(pd.DataFrame([[exp_ale_tuple[0],
                                                                        exp_ale_tuple[1],
                                                                        mut_type,
                                                                        mut_type_fraction]],
                                                              columns=["exp",
                                                                       "ale",
                                                                       df_column_name,
                                                                       "fraction"]))
    return mut_type_frac_df


def get_exp_avg_ale_mut_type_frac_df(mut_df, column_name):
    ale_mut_type_frac_df = get_ALE_mut_type_frac_df(mut_df, column_name)
    exp_avg_ale_mut_type_frac_df = pd.DataFrame()
    for exp_name_mut_type_tuple, mut_type_df in ale_mut_type_frac_df.groupby(["exp", column_name]):
        exp_name = exp_name_mut_type_tuple[0]
        mut_type = exp_name_mut_type_tuple[1]
        mut_type_ale_mean = np.mean(mut_type_df["fraction"])
        df = pd.DataFrame([[exp_name, mut_type, mut_type_ale_mean]],
                          columns=["experiment", column_name, "fraction"])
        exp_avg_ale_mut_type_frac_df = exp_avg_ale_mut_type_frac_df.append(df)
    return exp_avg_ale_mut_type_frac_df    
'''


def get_mut_type_frac_across_class_df(mut_df, exp_level_l, mut_type_class, all_mut_type_class_set):
    dfs = []  # Create a list to store DataFrames
    for class_tup, exp_ale_max_mut_df in mut_df.groupby(exp_level_l):
        for mut_type in all_mut_type_class_set:
            mut_type_count = len(exp_ale_max_mut_df[exp_ale_max_mut_df[mut_type_class]==mut_type])
            mut_type_frac = mut_type_count/len(exp_ale_max_mut_df)
            
            if len(exp_level_l) > 1:
                l = list(class_tup)
            else:
                l = [class_tup]
            data_l = l + [mut_type, mut_type_frac]
            col_name_l = exp_level_l + [mut_type_class, "fraction"]

            dfs.append(pd.DataFrame([data_l], columns=col_name_l))
    
    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()


def get_mut_type_avg_frac_across_class_df(mut_df, exp_level_l, mut_type_class, all_mut_type_class_set):
    mut_type_frac_df = get_mut_type_frac_across_class_df(mut_df, exp_level_l, mut_type_class, all_mut_type_class_set)
    dfs = []  # Create a list to store DataFrames
    
    # !!! Assuming the first class is the final class to consider penetration across
    final_penetration_class = exp_level_l[0]
    for tup, mut_type_df in mut_type_frac_df.groupby([final_penetration_class, mut_type_class]):
        final_class_type = tup[0]
        mut_type = tup[1]
        mut_type_class_type_mean = np.mean(mut_type_df["fraction"])
        df = pd.DataFrame([[final_class_type, mut_type, mut_type_class_type_mean]],
                         columns=[final_penetration_class, mut_type_class, "fraction"])
        dfs.append(df)
    
    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()
