import random
import string
import pandas as pd
from enum import Enum, auto
from itertools import chain
from typing import Optional
from biocypher._logger import logger

logger.debug(f"Loading module {__name__}.")


class WheatomicsAdapterNodeType(Enum):
    """
    Define types of nodes the adapter can provide.
    """

    GENE = "Gene"
    TRANSCRIPT = "Transcript"
    ONTOLOGY_TERM = "GO_term"
    # OFFICER = ":Officer"
    # LOCATION = ":Location"
    # CRIME = ":Crime"
    # PHONE_CALL = ":PhoneCall"
    # OBJECT = ":Object"


class WheatomicsAdapterGeneField(Enum):
    """
    Define possible fields the adapter can provide for genes.
    """
# MERGE (g:Gene {identifier: row.gene_id, chromosome: row.chrom, start: row.start, stop: row.stop, annotation: "Renan_refseqv2.1", genotype: "Renan" })

    ID = "_id"
    IDENTIFIER = "identifier"
    CHROMOSOME = "chromosome"
    START = "start"
    STOP = "stop"
    ANNOTATION = "annotation"
    GENOTYPE = "genotype"
       
class WheatomicsAdapterOntologyTermField(Enum):
    """
    Define possible fields the adapter can provide for annotations.
     """
    ID = "_id"
    GO_TERM = "GO_term"
    DESCRIPTION = "description"
    LABELS = "GO_term"
    
class WheatomicsAdapterTranscriptField(Enum): #TODO update
    """
    Define possible fields the adapter can provide for officers.
    """

# Transcripts
# - position
# - abundance value
# - gene-transcript-rel
# - condition
# - metadata.tsv
# 	- development_stage
# 	- Organ
# 	- Traitement (stress)

    ID = "_id"
    ABUNDANCE = "abundance"
# TODO complete

# Annotation
# - gene_id
# - 
# class WheatomicsAdapterOfficerField(Enum): #TODO update
#     """
#     Define possible fields the adapter can provide for officers.
#     """

#     ID = "_id"
#     NAME = "name"
#     SURNAME = "surname"
#     RANK = "rank"
#     BADGE_NUMBER = "badge_no"


class WheatomicsAdapterEdgeType(Enum):
    """
    Define possible edges the adapter can provide.
    """

    # HAS_TRANSCRIPT = "HAS_TRANSCRIPT" # BETWEEN GENE AND TRANSCRIPT
    HOMOLOGOUS_TO = "HOMOLOGOUS_TO" # BETWEEN 2 GENES
    LABELS = "HOMOLOGOUS_TO"
    RELATED_TO = "RELATED_TO" # BETWEEN GENE AND ANNOTATION

    # MAPPED = "MAPPED" # BETWEEN cs GENES AND GENES OF OTHER ACCESSIONS
    # KNOWS = "KNOWS"
    # INVOLVED_IN = "INVOLVED_IN"
    # MADE_CALL = "CALLED"  # abstract
    # RECEIVED_CALL = "CALLER"  # abstract
    # OCCURED_AT = "OCCURRED_AT"
    # INVESTIGATED_BY = "INVESTIGATED_BY"
    # PARTY_TO = "PARTY_TO"
    # RELATED_TO = "FAMILY_REL"


class WheatomicsAdapter:
    """
    Example BioCypher adapter. Generates nodes and edges for creating a
    knowledge graph.

    Args:
        node_types: List of node types to include in the result.
        node_fields: List of node fields to include in the result.
        edge_types: List of edge types to include in the result.
        edge_fields: List of edge fields to include in the result.
    """

    def __init__(
        self,
    ):
        __init__(self,None,None,None,None)
    

    def __init__(
        self,
        node_types: Optional[list] = None,
        node_fields: Optional[list] = None,
        edge_types: Optional[list] = None,
        edge_fields: Optional[list] = None,
    ):
        self._set_types_and_fields(node_types, node_fields, edge_types, edge_fields)
        
        # Define a constant for the columns to be renamed
        RENAME_COLUMNS = {
            "Gene stable ID": "gene_id",
            "Chromosome/scaffold name": "chromosome",
            "Gene start (bp)": "start",
            "Gene end (bp)": "stop"
        }

        rice_genes = self._read_gene_csv(csv_file="data/gene/mart_osjaponica_irgsp1_export.txt")
        rice_genes.drop(['Gene description','Gene type','UniProtKB Gene Name ID','UniProtKB Gene Name symbol'], axis=1, inplace=True)
        rice_genes.rename(columns=RENAME_COLUMNS, inplace=True)
        rice_genes['_id'] = rice_genes['gene_id']
        rice_genes = rice_genes.assign(annotation="irgsp1", genotype="Oryza sativa japonica")

        ath_genes = self._read_gene_csv(csv_file="data/gene/mart_ath_tair10_export.txt")
        ath_genes.drop(['Gene name', 'Gene type'], axis=1, inplace=True)
        ath_genes.rename(columns=RENAME_COLUMNS, inplace=True)
        ath_genes['_id'] = ath_genes['gene_id']
        ath_genes = ath_genes.assign(annotation="tair10", genotype="Col-0")

        csv1_wheat_genes = self._read_gene_csv(csv_file="data/gene/mart_tritaestivum_CS_IWGSCv1.1_export.txt")
        csv1_wheat_genes.drop(labels=['Gene name','Gene description'],axis=1, inplace=True)
        csv1_wheat_genes.rename(columns=RENAME_COLUMNS, inplace=True)
        csv1_wheat_genes['_id'] = csv1_wheat_genes['gene_id']
        csv1_wheat_genes = csv1_wheat_genes.assign(annotation="CSv1.1", genotype="Chinese Spring")

        renan_wheat_genes = self._read_gene_csv(csv_file="data/gene/TaeRenan_refseqv2.1_genesHC.tsv")
        renan_wheat_genes.rename(columns={"chrom": "chromosome"}, inplace=True)
        renan_wheat_genes['_id'] = renan_wheat_genes['gene_id']
        renan_wheat_genes = renan_wheat_genes.assign(annotation="Renan", genotype="Renan")

        self._node_data = pd.concat([rice_genes, ath_genes, csv1_wheat_genes, renan_wheat_genes])
        
        self._go_terms = self._get_annotation_data()
        self._edge_data = self._read_homolog_csv()
        self._oryzabase_annotations = self._read_oryzabase_csv()
        #print(self._node_data.head())
        #print(self._node_data.tail())
        # print head of the go_terms data frame
        print(self._go_terms.head())
        self._goa = self._read_goa_csv()
        print(self._goa.head())
        # print head of the edge data frame
        #print(self._edge_data.head())
        # self._data_homologs = self._read_csv(csv_file="homology.csv")
        # self._data = self._read_csv()
        # self._node_data = self._get_node_data()
        # self._edge_data = self._get_edge_data()
        
        # # print unique _labels
        # print(f"Unique labels: {self._data['_labels'].unique()}")

        # # print unique _type
        # print(f"Unique types: {self._data['_type'].unique()}")

    def _read_homolog_csv(self, csv_file='data/homology/Wheat_othologs_with_arabido_and_O.Sativa.japonic.txt'):
        """
        Read data from CSV file.
        """

        logger.info("Reading homolog data from CSV file.")
        df = pd.read_csv(csv_file, dtype=str, sep='\t')

        # Initialize a list to store the tuples (+ the relationship type)
        gene_pairs = []

        # Iterate over the DataFrame rows
        for index, row in df.iterrows():
            # Check if 'Ath gene' or 'Rice gene' exists and add the corresponding pair
            if pd.notna(row['Arabidopsis thaliana gene stable ID']):
                gene_pairs.append((row['Gene stable ID'], row['Arabidopsis thaliana gene stable ID'], "HOMOLOGOUS_TO"))
            if pd.notna(row['Oryza sativa Japonica Group gene stable ID']):
                gene_pairs.append((row['Gene stable ID'], row['Oryza sativa Japonica Group gene stable ID'], "HOMOLOGOUS_TO"))

        # Create a new DataFrame from the gene pairs
        result_df = pd.DataFrame(gene_pairs, columns=['source', 'target', '_type'])
        return result_df
 
    def _get_annotation_data(self):
        """
        Get GO terms from _read_annotation()
        """
        logger.info("Reading annotation data.")
        # get the data from _read_annotation()
        # as a pandas data frame
        data = self._read_annotation()
        # push result to a dictionary
        go_terms = []
        # loop through the data frame and create a new data frame
        for index, row in data.iterrows():
            # get the GO terms from the row
            GO_string = str(row["Gene Ontology"])
            if GO_string != 'nan':
                GO_terms = self._parse_go_terms(GO_string)
                for GO_term in GO_terms:
                    go_term = self._parse_go_term(GO_term)
                    DESC = self._parse_go_term_description(GO_term)
                    # create a new data frame with the GO_term, and description
                    go_terms.append({'_id': go_term,'GO_term': go_term, 'description': DESC, '_labels':'GO_term'})
                
        return pd.DataFrame(go_terms)
    
 
    def _read_annotation(self):
        """
        Read data from CSV file.
        """
        logger.info("Reading annotation data.")
        data = pd.read_csv('data/goa/OryzabaseGeneListEn_20241017010108.txt', sep='\t', delimiter=None, dtype='str', skip_blank_lines=True)
        data.drop(['CGSNL Gene Symbol','Gene symbol synonym(s)','CGSNL Gene Name','Gene name synonym(s)',\
            'Protein Name','Allele','Chromosome No.','Explanation','Trait Class','Gramene ID','Arm','Locate(cM)'], axis=1, inplace=True) 
        return data
    
    def _get_homolog_data(self):
        logger.info("Reading homolog data.")
        self._data = self._read_homolog_csv()
        
    def _read_oryzabase_csv(self):
        """
        Read data from CSV file.
        """
        logger.info(f"Reading oryzabase annotation data from CSV file.")
       # call _read_annotation() to get the data
        data = self._read_annotation()
        # data = pd.read_csv(csv_file, sep='\t', delimiter=None, dtype='str', skip_blank_lines=True)
        # data.fillna('', inplace=True)
        # data.drop(['CGSNL Gene Symbol','Gene symbol synonym(s)','CGSNL Gene Name','Gene name synonym(s)',\
            # 'Protein Name','Allele','Chromosome No.','Explanation','Trait Class','Gramene ID','Arm','Locate(cM)'], axis=1, inplace=True)
        # loop through the data frame and yield each row as a triple of gene_id, field, value
        for index, row in data.iterrows():
            gene_id = str(row["RAP ID"])
            go_annotations = []
            if not gene_id:
                continue
            else:
                GO_string = str(row["Gene Ontology"])
                if GO_string != 'nan':
                    GO_terms = self._parse_go_terms(GO_string)
                    for GO_term in GO_terms:
                        GO_term = self._parse_go_term(GO_term)
                        go_annotations.append({'source': gene_id, 'target': GO_term, '_type': 'RELATED_TO'})
                else:
                    continue
        return pd.DataFrame(go_annotations)
    
    def _parse_go_terms(self, go_terms):
        """
        Parse GO terms from a string.
        """
        return go_terms.split(';')
    
    def _parse_go_term(self, go_term):
        """
        Parse GO term from a string.
         """
        go_term = go_term.split('-')[0]
        return go_term.replace(' ', '')
   # Parse GO term description from a string
    def _parse_go_term_description(self, go_term):
        """
      Parse GO term description from a string.
         """
        description = go_term.split('-')[-1].strip()
        description = description.replace("'", "") 
        return description.replace('   ', '')
    
      # read the annotation data goa file
    # and return a pandas dataframe
    def _read_goa_csv(self, csv_file='data/goa/tair_annotations.gaf'):
        """
        Read data from CSV file.
        """
        logger.info("Reading annotation data.")
        # read the csv file
        # skip line starting with '!'
        # define names of the columns automatically with a prefix
        # define the data types of each column
        # skip blank lines  
        # skip lines starting with '!'
        goa = pd.read_csv(csv_file, sep='\t', delimiter=None, dtype='str', skip_blank_lines=True, comment='!', header=None)
        goa.rename('col_{}'.format, axis=1, inplace=True)
        # remove all columns except col_1 and col_4
        goa_clean = goa.loc[:, goa.columns.intersection(['col_1','col_4'])]
        goa_clean = goa_clean.assign(_labels='GO_term') # add new _labels column with Gene value 
        goa_clean['GO_term'] = goa_clean['col_4']
        # rename the columns col_4 to _id
        goa_clean.rename(columns={"col_4": "_id"}, inplace=True)
        # add a new column description with the value empty string
        goa_clean['description'] = ''
        return goa_clean
    
    def _read_gene_csv(self, csv_file):
        """
        Read data from CSV file.
        """
#         Gene stable ID	Chromosome/scaffold name	Gene start (bp)	Gene end (bp)	Gene type	Gene description	UniProtKB Gene Name symbol	UniProtKB Gene Name ID
# Os02g0788600	2	33502431	33505888	protein_coding			
        logger.info(f"Reading gene data from CSV file: {csv_file}.")

        data = pd.read_csv(csv_file, dtype=str, sep='\t')
        data = data.assign(_labels='Gene') # add new _labels column with Gene value

        # screen the entire data frame for double quotes
        data = data.map(lambda x: x.replace('"', "") if isinstance(x, str) else x)

        # rename 'id' to 'hash'
        #data.rename(columns={"Gene stable ID": "hash"}, inplace=True)

        return data

    def _get_node_data(self):
        """
        Get all rows that do not have a _type.
        """
        return self._data[self._data["_type"].isnull()]

    def _get_edge_data(self):
        """
        Get all rows that have a _type.
        """
        return self._data[self._data["_type"].notnull()]


    def _get_caller_data(self):
        """
        Subset to only call initiator relationships.
        """
        return self._data[self._data["_type"] == "CALLER"][["_start", "_end"]]

    def _get_called_data(self):
        """
        Subset to only call receiver relationships.
        """
        return self._data[self._data["_type"] == "CALLED"][["_start", "_end"]]

    def _get_phone(self, _id):
        """
        Get phone number for person.
        """
        if not _id in self._phone_data["_start"].values:
            return None

        phone_id = self._phone_data[self._phone_data["_start"] == _id]["_end"].values[0]
        phone = self._data[self._data["_id"] == phone_id]["phoneNo"].values[0]
        return phone

    def _get_email(self, _id):
        """
        Get email address for person.
        """
        if not _id in self._email_data["_start"].values:
            return None

        email_id = self._email_data[self._email_data["_start"] == _id]["_end"].values[0]
        email = self._data[self._data["_id"] == email_id]["email_address"].values[0]
        return email

    def get_nodes(self):
        """
        Returns a generator of node tuples for node types specified in the
        adapter constructor.
        """

        logger.info("Generating nodes.")
        # nodes: tuples of id, type, fields
        it=self._node_data.iterrows()
        for index, row in it:
            if row["_labels"] not in self.node_types:
                continue

            _id = row["_id"]
            _type = row["_labels"]
            _props = row.to_dict()
            # could filter non-values here

            # special cases - processing
            # if _type == WheatomicsAdapterNodeType.GENE.value:
            #     _props["phone"] = self._get_phone(_id)
            #     _props["email"] = self._get_email(_id)
            yield (
                _id,
                _type,
                _props,
            )
        for index, row in self._go_terms.iterrows():
            if row["_labels"] not in self.node_types:
                continue
            
            _id = row["_id"]
            _type = row["_labels"]
            _props = row.to_dict()
            yield (
                _id,
                _type,
                _props,
            )
        for index, row in self._goa.iterrows():
            if row["_labels"] not in self.node_types:
                continue
            
            _id = row["_id"]
            _type = row["_labels"]
            _props = row.to_dict()
            yield (
                _id,
                 _type,
                 _props,
             )

    def get_edges(self):
        """
        Returns a generator of edge tuples for edge types specified in the
        adapter constructor.

        Args:
            probability: Probability of generating an edge between two nodes.
        """

        logger.info("Generating edges.")
        for index, row in self._edge_data.iterrows():
            # if row["_type"] not in self.edge_types:
            #     continue
            _id = None
            _props = {}
            yield(
                _id,
                row['source'],
                row['target'],
                row['_type'],
                _props
            )
        for index, row in self._oryzabase_annotations.iterrows():
            # if row["_type"] not in self.edge_types:
            #     continue
            _id = None
            _props = {}
            yield(
                _id,
                row['source'],
                row['target'],
                row['_type'],
                _props
            )

    def _set_types_and_fields(self, node_types, node_fields, edge_types, edge_fields):
        if node_types:
            self.node_types = [type.value for type in node_types]
        else:
            self.node_types = [type.value for type in WheatomicsAdapterNodeType]

        if node_fields:
            self.node_fields = [field.value for field in node_fields]
        else:
            self.node_fields = [
                field.value
                for field in chain(
                    WheatomicsAdapterGeneField,
                    # WheatomicsAdapterTranscriptField,
                )
            ]

        if edge_types:
            self.edge_types = [type.value for type in edge_types]
        else:
            self.edge_types = [type.value for type in WheatomicsAdapterEdgeType]

        if edge_fields:
            self.edge_fields = [field.value for field in edge_fields]
        else:
            self.edge_fields = [field.value for field in chain()]
