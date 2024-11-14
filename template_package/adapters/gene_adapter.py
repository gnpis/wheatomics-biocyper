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

    GENE = ":Gene"
    TRANSCRIPT = ":Transcript"
    LABELS = "Gene"
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
        rice_genes = self._read_gene_csv(csv_file="data/gene/mart_osjaponica_irgsp1_export.txt")
        rice_genes.drop(['Gene description','Gene type','UniProtKB Gene Name ID','UniProtKB Gene Name symbol'], axis=1, inplace=True)
        ath_genes = self._read_gene_csv(csv_file="data/gene/mart_ath_tair10_export.txt")
        ath_genes.drop(['Gene name', 'Gene type'],axis=1, inplace=True)
        wheat_genes = self._read_gene_csv(csv_file="data/gene/mart_tritaestivum_CS_IWGSCv1.1_export.txt")
        wheat_genes.drop(labels=['Gene name','Gene description'],axis=1, inplace=True)
        self._node_data = pd.concat([rice_genes, ath_genes,wheat_genes])

        self._edge_data = self._read_homolog_csv()
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
        data.rename(columns={"id": "hash"}, inplace=True)

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
        iter=self._node_data.iterrows()
        for index, row in iter:
            if row["_labels"] not in self.node_types:
                continue

            _id = row["Gene stable ID"]
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

    def get_edges(self):
        """
        Returns a generator of edge tuples for edge types specified in the
        adapter constructor.

        Args:
            probability: Probability of generating an edge between two nodes.
        """

        logger.info("Generating edges.")
        for row in self._edge_data:
            # if row["_type"] not in self.edge_types:
            #     continue
            _id = None
            _props = {}
            yield(
                _id,
                row[0],
                row[1],
                row[2],
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
