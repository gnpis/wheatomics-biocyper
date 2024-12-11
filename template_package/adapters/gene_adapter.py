import random
import os
import re
import string
import time
import glob
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
    LABELS = "Gene"
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

    ID = "_id"
    COVERAGE = 'Coverage'
    FPKM = "FPKM"
    REFERENCE = "Reference"
    START = "Start"
    END = "End"
    STRAND = "Strand"
    TPM = 'TPM'
    ANNOTATION = "annotation"
    GENOTYPE = "genotype"
    # TODO: add missing properties for rnaseq:
    # - condition
    # - metadata.tsv
    # 	- development_stage
    # 	- Organ
    # 	- Traitement (stress)


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

    HAS_TRANSCRIPT = "HAS_TRANSCRIPT" # BETWEEN GENE AND TRANSCRIPT
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

        rice_genes = self.process_gene_data(
            csv_file="data/gene/mart_osjaponica_irgsp1_export.txt",
            drop_columns=['Gene description', 'Gene type', 'UniProtKB Gene Name ID', 'UniProtKB Gene Name symbol'],
            rename_columns=RENAME_COLUMNS,
            annotation="irgsp1",
            genotype="Oryza sativa japonica"
        )

        ath_genes = self.process_gene_data(
            csv_file="data/gene/mart_ath_tair10_export.txt",
            drop_columns=['Gene name', 'Gene type'],
            rename_columns=RENAME_COLUMNS,
            annotation="tair10",
            genotype="Col-0"
        )

        csv1_wheat_genes = self.process_gene_data(
            csv_file="data/gene/mart_tritaestivum_CS_IWGSCv1.1_export.txt",
            drop_columns=['Gene name', 'Gene description'],
            rename_columns=RENAME_COLUMNS,
            annotation="CSv1.1",
            genotype="Chinese Spring"
        )

        # csv2_wheat_genes = pd.DataFrame()
        csv2_wheat_genes = self.process_wheat_csv2_gene_data(
            csv_file="data/gene/iwgsc_refseqv2.1_annotation_200916_HC.gff3",
            annotation="CSv2.1",
            genotype="Chinese spring"
        )

        renan_wheat_genes = self.process_gene_data(
            csv_file="data/gene/TaeRenan_refseqv2.1_genesHC.tsv",
            drop_columns=[],
            rename_columns={"chrom": "chromosome"},
            annotation="Renan",
            genotype="Renan"
        )

        # TODO: fix some glitches in rice homologs:
        # TraesCS7A02G230600		gene-rps8	HOMOLOGOUS_TO

        self._node_genes_data = pd.concat([rice_genes, ath_genes, csv1_wheat_genes, renan_wheat_genes, csv2_wheat_genes])

        self._edge_data = self._read_homolog_csv()
        wheat_homologs_df = self._read_wheat_homolog_csv()
        self._edge_data = pd.concat([self._edge_data, wheat_homologs_df])

        start_time = time.time()
        self._node_data_rnaseq = pd.DataFrame()
        self._edge_gene_transcripts_df = pd.DataFrame()

        gene_id_column = 'Gene ID'
        # rnaseq_df=self._read_transcript_csv('data/transcript/')
        for rnaseq_df in self._read_transcript_csv('data/transcript/'):
            logger.info(f"rnaseq_df size is: {len(rnaseq_df)}")
            columns_to_drop = ['Start', 'End', 'FPKM','Strand', 'TPM']
            rnaseq_df_dropped = rnaseq_df.drop(columns=columns_to_drop) # drop useless columns for node, this data is used only for the relationship/edge
            self._node_data_rnaseq = pd.concat([rnaseq_df_dropped, self._node_data_rnaseq])
            # rnaseq_df=rnaseq_df.assign(_type=WheatomicsAdapterEdgeType.HAS_TRANSCRIPT.value)
            rnaseq_df.rename(columns={'ref_gene_id':'source',gene_id_column:'target'}, inplace=True)
            self._edge_gene_transcripts_df = pd.concat([rnaseq_df, self._edge_gene_transcripts_df])
            logger.info(f"TEMP self._edge_gene_transcripts_df size is: {len(self._edge_gene_transcripts_df)}")
        self._edge_gene_transcripts_df = self._edge_gene_transcripts_df.assign(_type=WheatomicsAdapterEdgeType.HAS_TRANSCRIPT.value)

        logger.info(f"ALL self._edge_gene_transcripts_df size is: {len(self._edge_gene_transcripts_df)}")

        self._node_data_rnaseq = self._node_data_rnaseq.drop_duplicates(subset=[gene_id_column]).sort_values(by=[gene_id_column])
        self._node_data_rnaseq = self._node_data_rnaseq.assign(_labels=WheatomicsAdapterNodeType.TRANSCRIPT.value)

        end_time = time.time()
        execution_time = end_time - start_time
        logger.info(f"The _read_transcript_csv method took {execution_time} seconds to execute.")
        
        self._go_terms = self._get_annotation_data()
        self._goa_go_terms = self._get_goa_go_term()
        
        
        self._oryzabase_annotations = self._read_oryzabase_csv()
        self._annotation_edges = self._get_goa_annotation()

    def process_wheat_csv2_gene_data(self, csv_file, genotype, annotation):

        gff_df = pd.read_csv(
            csv_file,
            sep="\t",
            comment="#",
            header=None,
            names=["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
        )

        # Filter only rows where the type is 'gene'
        genes_df = gff_df[gff_df['type'] == 'gene'].copy()
        
        # Extract relevant information from the attributes column for the ID
        # genes_df['ID'] = genes_df['attributes'].str.extract(r'ID=([^;]+)')
        genes_df.loc[:, 'ID'] = genes_df['attributes'].str.extract(r'ID=([^;]+)')
        genes_df['_id'] = genes_df['ID']
        
        # gene_id chromosome      start   stop    _labels _id     annotation      genotype        ID      Chromosome
        # Os02g0788600    2       33502431        33505888        Gene    Os02g0788600    irgsp1  Oryza sativa japonica

        # Create the desired output DataFrame
        output_df = genes_df[['ID', 'seqid', 'start', 'end', '_id']].rename(columns={
            "ID": "gene_id",
            "seqid": "chromosome",
            "start": "start",
            "end": "stop"
        }).assign(annotation=annotation,genotype=genotype,_labels=WheatomicsAdapterNodeType.GENE.value) # add new _labels column with Gene value

        return output_df

    def process_gene_data(self, csv_file, drop_columns, rename_columns, annotation, genotype):
        genes = self._read_gene_csv(csv_file=csv_file)
        if drop_columns:
            genes.drop(drop_columns, axis=1, inplace=True)
        genes.rename(columns=rename_columns, inplace=True)
        genes['_id'] = genes['gene_id']
        genes = genes.assign(annotation=annotation, genotype=genotype)
        return genes

    def _find_file_pairs(self, directory):
        # Get a list of all files ending with 'refmap' or 'abundance.tsv' in the given directory
        files = [os.path.join(directory, f) for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f)) and (f.endswith("refmap") or f.endswith("abundance.tsv"))]
        sorted_files=sorted(files)
        # Find pairs of files with the same base name
        file_pairs = []
        for file in sorted_files:
            file_pairs.append(file)
            if len(file_pairs) == 2:
                yield(file_pairs)
                file_pairs = []
                # break # TODO remove, for debug only

    def _read_transcript_csv(self, csv_dir):
        """
        Read transcript data from CSV files.
        """

        logger.info(f"Reading data from CSV files in {csv_dir}.")
        refmap_dict = {}

        # Function to replace Gene ID with ref_id using the mapping
        def map_gene_id_to_ref_info(gene_id):
            return refmap_dict.get(gene_id, (gene_id, None))  # If not found, return the original gene_id and None for ref_gene_id

        file_pairs = self._find_file_pairs(csv_dir)
        for pair in file_pairs:
            abundance_df = pd.read_csv(pair[0], sep='\t')
            refmap_df = pd.read_csv(pair[1], sep='\t')

            # Create a mapping of the Gene ID to ref_id from refmap
            for _, row in refmap_df.iterrows():
                qry_ids = row['qry_id_list'].split(',')  # Split the qry_id_list by commas
                for qry_id in qry_ids:
                    # Take the part before the first pipe as the key (Gene ID) and the ref_id as the value
                    gene_id = qry_id.split('|')[0]
                    refmap_dict[gene_id] = (row['ref_id'], row['ref_gene_id'])

            # Apply the mapping function to the 'Gene ID' column in the abundance DataFrame
            abundance_df[['Gene ID', 'ref_gene_id']] = abundance_df['Gene ID'].apply(
                lambda gene_id: pd.Series(map_gene_id_to_ref_info(gene_id))
            )
            abundance_df = abundance_df.dropna(subset=['ref_gene_id'])
            abundance_df=abundance_df.assign(_labels=WheatomicsAdapterNodeType.TRANSCRIPT.value)
            abundance_df["_id"]=abundance_df["Gene ID"]

            yield abundance_df

    def _read_csv(self, csv_file):
        """
        Read data from CSV file.
        """

        logger.info("Reading homolog data from CSV file.")

        data = pd.read_csv(csv_file, dtype=str, sep='\t')

        for index, row in data.iterrows():
            # if row["_labels"] not in self.node_types:
            #     continue

            wheat_gene_stable_id = str(row["Gene stable ID"])
            arabidopsis_gene_id = str(row["Arabidopsis thaliana gene stable ID"])
            japonica_gene_id = str(row["Oryza sativa Japonica Group gene stable ID"])

            if arabidopsis_gene_id != 'nan':
                yield (
                    wheat_gene_stable_id,
                    arabidopsis_gene_id,
                    'HOMOLOGOUS_TO'
                )
            if japonica_gene_id != 'nan':
                yield (
                    wheat_gene_stable_id,
                    japonica_gene_id,
                    'HOMOLOGOUS_TO'
                )

        # rename 'id' to 'hash'
        data.rename(columns={"id": "hash"}, inplace=True)

        return data

    def _read_homolog_csv(self, csv_file='data/homology/Wheat_othologs_with_arabido_and_O.Sativa.japonic.txt'):
        """
        Read data from CSV file.
        """

        logger.info(f"Reading homolog data from CSV file: {csv_file}.")
        df = pd.read_csv(csv_file, dtype=str, sep='\t')

        # Start the timer
        start_time = time.time()
        
        # Define constants
        ATH_GENE_ID = 'Arabidopsis thaliana gene stable ID'
        OS_GENE_ID = 'Oryza sativa Japonica Group gene stable ID'
        WHEAT_GENE_ID = 'Gene stable ID'

        # Create masks for the 'Arabidopsis thaliana gene stable ID' and 'Oryza sativa Japonica Group gene stable ID' columns
        mask_ath = df[ATH_GENE_ID].notna()
        mask_os = df[OS_GENE_ID].notna()

        # Create a new DataFrame from the gene pairs using the mask
        result_df = pd.concat([
            df.loc[mask_ath, [WHEAT_GENE_ID, ATH_GENE_ID]].rename(columns={ATH_GENE_ID: 'target'}).assign(_type=WheatomicsAdapterEdgeType.HOMOLOGOUS_TO.value),
            df.loc[mask_os, [WHEAT_GENE_ID, OS_GENE_ID]].rename(columns={OS_GENE_ID: 'target'}).assign(_type=WheatomicsAdapterEdgeType.HOMOLOGOUS_TO.value)
        ]).rename(columns={WHEAT_GENE_ID: 'source'}).reset_index(drop=True)

        # Calculate the execution time
        end_time = time.time()
        execution_time = end_time - start_time
        logger.info(f"The method took {execution_time} seconds to execute.")

        return result_df
    
    def _read_wheat_homolog_csv(self, csv_file='data/homology/TaeRenan_refseqv2.1_CORRESPONDANCE_CSv1_CSv2.txt'):
        """
        Read data from CSV file.
        """

        logger.info(f"Reading homolog data from CSV file: {csv_file}.")
        df = pd.read_csv(csv_file, dtype=str, sep='\t')

        # Start the timer
        start_time = time.time()

        # Create a new DataFrame with the expanded rows
        result_df = pd.DataFrame()

        # First part: Keep csv2_1 and add the renan (csv1_1)
        result_df = pd.concat([result_df, df[['csv2_1', 'renan']].assign(_type=WheatomicsAdapterEdgeType.HOMOLOGOUS_TO.value).rename(columns={'renan': 'target'})])

        # Second part: Keep csv2_1 and add the csv1_1 value
        result_df = pd.concat([result_df, df[['csv2_1', 'csv1_1']].assign(_type=WheatomicsAdapterEdgeType.HOMOLOGOUS_TO.value).rename(columns={'csv1_1': 'target'})])

        # Rename the 'csv2_1' column to 'source'
        result_df.rename(columns={'csv2_1': 'source'}, inplace=True)
        end_time = time.time()
        execution_time = end_time - start_time
        logger.info(f"The method took {execution_time} seconds to execute.")

        return result_df

    def _read_gene_csv(self, csv_file):
        """
        Read data from CSV file.
        """
        #         Gene stable ID	Chromosome/scaffold name	Gene start (bp)	Gene end (bp)	Gene type	Gene description	UniProtKB Gene Name symbol	UniProtKB Gene Name ID
        # Os02g0788600	2	33502431	33505888	protein_coding			
        logger.info(f"Reading gene data from CSV file: {csv_file}.")

        data = pd.read_csv(csv_file, dtype=str, sep='\t')
        data = data.assign(_labels=WheatomicsAdapterNodeType.GENE.value) # add new _labels column with Gene value

        # screen the entire data frame for double quotes
        data = data.map(lambda x: x.replace('"', "") if isinstance(x, str) else x)

        # rename 'id' to 'hash'
        data.rename(columns={"id": "hash"}, inplace=True)

        return data
 
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
            go_string = str(row["Gene Ontology"])
            if go_string != 'nan':
                go_terms = self._parse_go_terms(go_string)
                for go_term in go_terms:
                    # create a new data frame with the GO_term, and description
                    go_terms.append({'_id': go_term[0],'GO_term': go_term[0], 'description': go_term[1].replace("'","prime"), '_labels':'GO_term'})

        return pd.DataFrame(go_terms)

    def _read_annotation(self):
        """
        Read data from CSV file.
        """
        oryza_annotation_file = "data/goa/OryzabaseGeneListEn_20241017010108.txt"
        logger.info(f"Reading annotation data from { oryza_annotation_file }.")
        data = pd.read_csv(oryza_annotation_file, sep='\t', delimiter=None, dtype='str', skip_blank_lines=True, escapechar='\\')
        data.drop(['CGSNL Gene Symbol','Gene symbol synonym(s)','CGSNL Gene Name','Gene name synonym(s)',\
            'Protein Name','Allele','Chromosome No.','Explanation','Trait Class','Gramene ID','Arm','Locate(cM)'], axis=1, inplace=True) 
        return data
        
    def _read_oryzabase_csv(self):
        """
        Read data from CSV file.
        """
        logger.info("Reading oryzabase annotation data from CSV file to create edges.")
        data = self._read_annotation()
        go_annotations = []
        for index, row in data.iterrows():
            gene_id = row["RAP ID"]
            if pd.isna(row["RAP ID"]) or pd.isnull(row["RAP ID"]):
                continue
            else:
                go_string = str(row["Gene Ontology"])
                if go_string != 'nan':
                    go_terms = self._parse_go_terms(go_string)
                    for go_term in go_terms:
                        go_annotations.append({'source': gene_id.strip(), 'target': go_term[0], '_type': 'RELATED_TO'})
                else:
                    continue
        
        df = pd.DataFrame(go_annotations)

        # Add double quotes to all string fields using apply with axis=1
        df = df.apply(lambda row: row.map(lambda x: f'"{x}"' if isinstance(x, str) else x), axis=1)
        return df.copy()
    
    def _parse_go_terms(self, go_terms):
        """
        Parse GO terms from a string using a Regular Expression:

        (GO:\d{7}):
            Matches the GO term (GO: followed by 7 digits).

        - (.*?):
            Matches the description part after the GO term and hyphen.
            The .*? ensures a non-greedy match to capture the description until the next GO term.

        (?=, GO:\d{7}|\Z):
            Positive lookahead to ensure the match stops at:
                , GO:\d{7}: The next GO term.
                \Z: The end of the string.
        """
        # Split the string into GO term and description pairs
        # pattern = r'(GO:\d{7}) - ([^,]+(?:, [^,]+)*)'
        pattern = r'(GO:\d{7}) - (.*?)(?=, GO:\d{7}|\Z)'
        
        matches = re.findall(pattern, go_terms)
        return matches
    
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
    
    def _read_goa_csv(self, csv_file='data/goa/tair_annotations.gaf'):
        """
        Read the annotation data goa file and return a pandas dataframe
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
        return goa_clean
    
    def _get_goa_go_term(self):
        # read the annotation data goa file
        # and return a pandas dataframe
        goa = self._read_goa_csv()
        goa = goa.assign(_labels='GO_term') # add new _labels column with Gene value 
        # drop col_1
        goa.drop(['col_1'], axis=1, inplace=True)
        goa['GO_term'] = goa['col_4']
        # rename the columns col_4 to _id
        goa.rename(columns={"col_4": "_id"}, inplace=True)
        # add a new column description with the value empty string
        goa['description'] = '' 
        return goa
    def _get_goa_annotation(self):
        annotation = self._read_goa_csv()
        # rename the columns col_1 to source and col_4 to target
        annotation.rename(columns={"col_1": "source", "col_4": "target"}, inplace=True)
        # add a new column _type with the value RELATED_TO
        annotation['_type'] = 'RELATED_TO'
        return annotation

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

        logger.info(f"Generating {len(self._node_genes_data)} gene nodes.")
        # nodes: tuples of id, type, fields
        
        gene_it=self._node_genes_data.iterrows()
        for index, row in gene_it:
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
        logger.info(f"Generating {len(self._node_data_rnaseq)} transcript nodes.")
        rnaseq_it=self._node_data_rnaseq.iterrows()
        for index, row in rnaseq_it:
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
        logger.info(f"Generating {len(self._go_terms)} go terms nodes.")
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
        logger.info(f"Generating {len(self._goa_go_terms)} goa go terms nodes.")
        for index, row in self._goa_go_terms.iterrows():
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

        logger.info(f"Generating {len(self._edge_data)} homology edges.")
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
        
        logger.info(f"Generating {len(self._edge_gene_transcripts_df)} rnaseq/transcript edges.")
        self._edge_gene_transcripts_df.to_csv('edge_gene_transcripts_df.csv', index=False)
        for index, row in self._edge_gene_transcripts_df.iterrows():
            # if row["_type"] not in self.edge_types:
            #     continue
            _props = row.to_dict()
            _id = str(hash((row['_id'], row['source'], row['target'], row['_type'], tuple(_props.items()))))
            if 'source' in row and 'target' in row: # ignore cases when transcripts are not linked to a structural gene (can be a ncRNA, or whatever found by the stringtie
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
        for index, row in self._annotation_edges.iterrows():
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
                    WheatomicsAdapterTranscriptField,
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
