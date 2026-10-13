/* Schema of a microbetag database built from your own genome catalogue, with microbetag-export-db.

   Load it in the MySQL container: place this file in /docker-entrypoint-initdb.d/ and the exported files in
   /var/lib/mysql-files/ (see server/services/db/README.md); MySQL runs it when it starts with an empty
   database. Same tables as init.sql (the GTDB database), with genome identifiers of any length, plus
   genomeTaxonomy, which the server uses to match taxa to the closest genomes.
*/

CREATE DATABASE IF NOT EXISTS microbetagDB;
USE microbetagDB;

-- Lineage of each genome of the catalogue
DROP TABLE IF EXISTS genomeTaxonomy;
CREATE TABLE genomeTaxonomy(
    genomeId VARCHAR(255),
    taxonomy TEXT,
    PRIMARY KEY (genomeId)
);
LOAD DATA INFILE '/var/lib/mysql-files/genome_taxonomy.tsv' INTO TABLE genomeTaxonomy
FIELDS TERMINATED BY '\t' LINES TERMINATED BY '\n';

-- phenotrex predictions: presence (YES/NO/N/A) and score of each trait
DROP TABLE IF EXISTS phenDB;
CREATE TABLE phenDB(
    gtdbId VARCHAR(255),
    aceticAcid VARCHAR(5),
    aceticAcidScore DECIMAL(5,4),
    aob VARCHAR(5),
    aobScore DECIMAL(5,4),
    aSaccharolytic VARCHAR(5),
    aSaccharolyticScore DECIMAL(5,4),
    autoCo2 VARCHAR(5),
    autoCo2Score DECIMAL(5,4),
    butanol VARCHAR(5),
    butanolScore DECIMAL(5,4),
    butyricAcid VARCHAR(5),
    butyricAcidScore DECIMAL(5,4),
    dGlucose VARCHAR(5),
    dGlucoseScore DECIMAL(5,4),
    dLacticAcid VARCHAR(5),
    dLacticAcidScore DECIMAL(5,4),
    ethanol VARCHAR(5),
    ethanolScore DECIMAL(5,4),
    fermentative VARCHAR(5),
    fermentativeScore DECIMAL(5,4),
    fixingN2 VARCHAR(5),
    fixingN2Score DECIMAL(5,4),
    formicAcid VARCHAR(5),
    formicAcidScore DECIMAL(5,4),
    halophilic VARCHAR(5),
    halophilicScore DECIMAL(5,4),
    hydrogen VARCHAR(5),
    hydrogenScore DECIMAL(5,4),
    indole VARCHAR(5),
    indoleScore DECIMAL(5,4),
    isobutyricAcid VARCHAR(5),
    isobutyricAcidScore DECIMAL(5,4),
    isovalericAcid VARCHAR(5),
    isovalericAcidScore DECIMAL(5,4),
    lLacticAcid VARCHAR(5),
    lLacticAcidScore DECIMAL(5,4),
    methanotroph VARCHAR(5),
    methanotrophScore DECIMAL(5,4),
    NOB VARCHAR(5),
    NOBScore DECIMAL(5,4),
    nonFermentative VARCHAR(5),
    nonFermentativeScore DECIMAL(5,4),
    phototrophy VARCHAR(5),
    phototrophyScore DECIMAL(5,4),
    psychrophilic VARCHAR(5),
    psychrophilicScore DECIMAL(5,4),
    rAcetoin VARCHAR(5),
    rAcetoinScore DECIMAL(5,4),
    saccharolytic VARCHAR(5),
    saccharolyticScore DECIMAL(5,4),
    succinicAcid VARCHAR(5),
    succinicAcidScore DECIMAL(5,4),
    sulfateReducer VARCHAR(5),
    sulfateReducerScore DECIMAL(5,4),
    symbiont VARCHAR(5),
    symbiontScore DECIMAL(5,4),
    T3SS VARCHAR(5),
    T3SSScore DECIMAL(5,4),
    T6SS VARCHAR(5),
    T6SSScore DECIMAL(5,4),
    thermophilic VARCHAR(5),
    thermophilicScore DECIMAL(5,4),
    anaerobe VARCHAR(5),
    anaerobeScore DECIMAL(5,4),
    aerobe VARCHAR(5),
    aerobeScore DECIMAL(5,4),
    PRIMARY KEY (gtdbId)
);
LOAD DATA INFILE '/var/lib/mysql-files/phen_predictions.tsv' INTO TABLE phenDB
FIELDS TERMINATED BY '\t' LINES TERMINATED BY '\n';

-- Pathway complementarities: each complement once, and the complements of each genome pair
DROP TABLE IF EXISTS uniqueComplements;
CREATE TABLE uniqueComplements(
    complementId INT,
    KoModuleId VARCHAR(10),
    complement TEXT,
    pathway TEXT,
    PRIMARY KEY (complementId)
);
LOAD DATA INFILE '/var/lib/mysql-files/unique_complements.tsv' INTO TABLE uniqueComplements
FIELDS TERMINATED BY '\t' LINES TERMINATED BY '\n';

DROP TABLE IF EXISTS pathwayComplementarity;
CREATE TABLE pathwayComplementarity(
    beneficiaryGenome VARCHAR(255),
    donorGenome VARCHAR(255),
    complmentId TEXT
);
LOAD DATA INFILE '/var/lib/mysql-files/pathway_complementarity.tsv' INTO TABLE pathwayComplementarity
FIELDS TERMINATED BY '\t' LINES TERMINATED BY '\n';
CREATE INDEX genome_pair ON pathwayComplementarity(beneficiaryGenome, donorGenome);

-- The GTDB database maps genomes to the PATRIC ids of their GEMs; here each genome is its own id
DROP TABLE IF EXISTS patricId2genomeId;
CREATE TABLE patricId2genomeId(
    gtdbGenomeAccession VARCHAR(255),
    patricId VARCHAR(255),
    PRIMARY KEY (patricId)
);
LOAD DATA INFILE '/var/lib/mysql-files/genome_ids.tsv' INTO TABLE patricId2genomeId
FIELDS TERMINATED BY '\t' LINES TERMINATED BY '\n';

-- Not used for custom catalogues (taxa are matched with genomeTaxonomy); kept empty for the server's queries
DROP TABLE IF EXISTS genome2taxNcbiId;
CREATE TABLE genome2taxNcbiId(
    ncbiTaxId VARCHAR(13),
    genomeId VARCHAR(255),
    PRIMARY KEY (genomeId)
);

-- Seed complementarities
DROP TABLE IF EXISTS seedScores;
CREATE TABLE seedScores(
    patricGenomeA VARCHAR(255),
    patricGenomeB VARCHAR(255),
    competitionScore DECIMAL(4,3),
    complementaritScore DECIMAL(4,3)
);
LOAD DATA INFILE '/var/lib/mysql-files/seed_scores.tsv' INTO TABLE seedScores
FIELDS TERMINATED BY '\t' LINES TERMINATED BY '\n';
CREATE INDEX genome_pair ON seedScores(patricGenomeA, patricGenomeB);

DROP TABLE IF EXISTS SeedComplements;
CREATE TABLE SeedComplements(
    patricBeneficaryId VARCHAR(255),
    patricDonorId VARCHAR(255),
    seedComplements TEXT,
    PRIMARY KEY (patricBeneficaryId, patricDonorId)
);
LOAD DATA INFILE '/var/lib/mysql-files/seed_complements.tsv' INTO TABLE SeedComplements
FIELDS TERMINATED BY '\t' LINES TERMINATED BY '\n';

DROP TABLE IF EXISTS nonSeedSets;
CREATE TABLE nonSeedSets(
    patricId VARCHAR(255),
    nonSeeds TEXT,
    PRIMARY KEY (patricId)
);
LOAD DATA INFILE '/var/lib/mysql-files/non_seed_sets.tsv' INTO TABLE nonSeedSets
FIELDS TERMINATED BY '\t' LINES TERMINATED BY '\n';
