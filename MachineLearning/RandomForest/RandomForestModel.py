
from MachineLearning.VEP_Annotation import Parser
from MachineLearning.VEP_Annotation import VEP_extract
import pandas as pd
import joblib
import os


def input(json_data):
    """
    :param json_data: The JSON data to be processed.
    :return: Processed DataFrame.
    """
    #print("Type of json_data:", type(json_data))
    df = VEP_extract.process_json_to_dataframe(json_data)
    # delete unwanted values
    df = df.drop(columns=['variant', 'transcript_id', 'clinvar_label'])
    df['lof_score'] = df['lof_score'].apply(Parser.parse_loftee)
    df['impact'] = df['impact'].apply(Parser.parse_impact)
    #List of classes for one-hot encoding
    classes = ['3_prime_UTR_variant', '5_prime_UTR_variant', 'coding_sequence_variant', 'downstream_gene_variant', 'frameshift_variant', 'inframe_deletion', 'inframe_insertion', 'intron_variant', 'mature_miRNA_variant', 'missense_variant', 'non_coding_transcript_exon_variant', 'protein_altering_variant', 'splice_acceptor_variant', 'splice_donor_5th_base_variant', 'splice_donor_region_variant', 'splice_donor_variant', 'splice_polypyrimidine_tract_variant', 'splice_region_variant', 'start_lost', 'stop_gained', 'stop_lost', 'stop_retained_variant', 'synonymous_variant', 'transcript_ablation', 'upstream_gene_variant']
    df_one_hot = Parser.single_hot_encoding(df['consequence'],classes)
    df = df.drop(columns=['consequence'])
    df = pd.concat([df, df_one_hot], axis=1)
    #print full frame in the console
    '''
    with pd.option_context('display.max_columns', None, 'display.width', None):
        print("DataFrame after processing:")
        print(df.head())
    
    '''
    array = df.to_numpy()
    return array 

# Import model from the file
from sklearn.ensemble import RandomForestRegressor
def predict(json_data):
    """
    Predicts the target variable using the Random Forest model.
    
    :param array: The input j_son data for prediction.
    :return: Predicted values in an array for example one sample as input,
             results in an array of length 1.
    """
    # Convert JSON data to array
    array = input(json_data)
    #regressor = joblib.load('MachineLearning/RandomForest/Random_forest_model.pkl')
    current_dir = os.path.dirname(os.path.abspath(__file__))
    model_path = os.path.join(current_dir, 'Random_forest_model.pkl')
    regressor = joblib.load(model_path)
    # Perform prediction using the loaded model
    predictions = regressor.predict(array)
    return predictions.tolist()