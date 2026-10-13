-- author: Haris Zafeiropoulos
-- email: haris.zafeiropoulos@kuleuven.be
-- description: Initiate microbetag database scheme and populate data
-- usage: mysql --local-infile=1 -u root -p < init.sql
-- notes: makes sure that you have first run SET GLOBAL local_infile=ON;
--        this will allow the copy of the .tsv/.csv files


SET GLOBAL local_infile=ON;
SET @@GLOBAL.local_infile = 1;
-- OPT_LOCAL_INFILE=1;

DROP DATABASE IF EXISTS microbetagDB;
CREATE DATABASE microbetagDB; 
USE microbetagDB;

/* ///////////////////////////////////////
	TABLE RELATED TO PHEND PREDICTIONS
/////////////////////////////////////// */

DROP TABLE IF EXISTS phenDB;
CREATE TABLE phenDB(
    gtdbId VARCHAR(15),
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
)
ENGINE=InnoDB 
ROW_FORMAT=COMPRESSED;

-- Disable keys and indexes
ALTER TABLE phenDB DISABLE KEYS;

LOAD DATA INFILE '/var/lib/mysql-files/gtdb_phen_predictions_oxy.tsv'
INTO TABLE phenDB
FIELDS TERMINATED BY '\t'
LINES TERMINATED BY '\n';

-- Enable again keys
ALTER TABLE phenDB ENABLE KEYS;



/* /////////////////////////////////////////////
	TABLEs RELATED TO PATHWAY COMPLEMENTARITIES
//////////////////////////////////////////////// */

-- Table for unique complementarity cases
DROP TABLE IF EXISTS uniqueComplements; 
CREATE TABLE uniqueComplements(
    complementId INT AUTO_INCREMENT,
    KoModuleId VARCHAR(7),
    complement VARCHAR(200),
    pathway VARCHAR(250),
    PRIMARY KEY (complementId)
);

LOAD DATA INFILE '/var/lib/mysql-files/unique_complements.tsv'
INTO TABLE uniqueComplements
FIELDS TERMINATED BY '\t'
LINES TERMINATED BY '\n';

-- remove new line characters
UPDATE uniqueComplements SET pathway = REPLACE(REPLACE(pathway, '\r', ''), '\n', '');


-- Table for GenomeId 2 NCBI TaxId 
DROP TABLE IF EXISTS genome2taxNcbiId; 
CREATE TABLE genome2taxNcbiId(
    ncbiTaxId VARCHAR(13),
	genomeId VARCHAR(20),
    PRIMARY KEY (genomeId)
);

LOAD DATA INFILE '/var/lib/mysql-files/genomeId2ncbiTaxId.tsv'
INTO TABLE genome2taxNcbiId
FIELDS TERMINATED BY '\t'
LINES TERMINATED BY '\n';


CREATE TABLE your_temporary_table(
    ncbiTaxId VARCHAR(13),
	genomeId VARCHAR(20),
    PRIMARY KEY (genomeId)
);

LOAD DATA INFILE '/var/lib/mysql-files/update_genomeId2ncbiTaxId.tsv' INTO TABLE your_temporary_table
FIELDS TERMINATED BY '\t' LINES TERMINATED BY '\n';

INSERT IGNORE INTO genome2taxNcbiId
SELECT * FROM your_temporary_table;

DROP TABLE your_temporary_table;

-- Table for pathway complementarities
DROP TABLE IF EXISTS pathwayComplementarity;
CREATE TABLE pathwayComplementarity(
	beneficiaryGenome VARCHAR(15),
	donorGenome VARCHAR(15),
	complmentId VARCHAR(1000)
);

LOAD DATA INFILE '/var/lib/mysql-files/complementarities/part_0_mapped_only_gtdb_genomes.tsv'
INTO TABLE pathwayComplementarity
FIELDS TERMINATED BY '\t'
LINES TERMINATED BY '\n';


CREATE INDEX genome_pair ON pathwayComplementarity(beneficiaryGenome, donorGenome);

-- LOAD DATA INFILE '/var/lib/mysql-files/complementarities/part_4_mapped_only_gtdb_genomes.tsv'
-- INTO TABLE pathwayComplementarity
-- FIELDS TERMINATED BY '\t'
-- LINES TERMINATED BY '\n'
-- IGNORE 1 LINES;



/* /////////////////////////////////////////////
	TABLEs RELATED TO SEED SCORES
//////////////////////////////////////////////// */
DROP TABLE IF EXISTS patricId2genomeId;
CREATE TABLE patricId2genomeId(
	gtdbGenomeAccession VARCHAR(15),
	patricId VARCHAR(12),
	PRIMARY KEY (patricId)
);


LOAD DATA INFILE '/var/lib/mysql-files/gtdb2patricIds.tsv'
INTO TABLE patricId2genomeId
FIELDS TERMINATED BY '\t'
LINES TERMINATED BY '\n';




DROP TABLE IF EXISTS seedScores;
CREATE TABLE seedScores(
	patricGenomeA VARCHAR(15),
	patricGenomeB VARCHAR(15),
	competitionScore DECIMAL(4,3),
	complementaritScore DECIMAL(4,3)
);


LOAD DATA INFILE '/var/lib/mysql-files/all_scores.tsv'
INTO TABLE seedScores
FIELDS TERMINATED BY '\t'
LINES TERMINATED BY '\n';


CREATE INDEX genome_pair ON seedScores(patricGenomeA, patricGenomeB);

/* remember!
run mysql like this: 
mysql --local-infile -u username -p
to enable LOCAL INFILE
*/



/* /////////////////////////////////////////////
	TABLE FOR SEED COMPLEMENTARITIES
//////////////////////////////////////////////// */



DROP TABLE IF EXISTS SeedComplements;
CREATE TABLE SeedComplements( 
	patricBeneficaryId VARCHAR(15), 
	patricDonorId VARCHAR(15), 
	seedComplements VARCHAR(1000), 
	PRIMARY KEY (patricBeneficaryId, patricDonorId)
);









CREATE TABLE nonSeedSets( 
	patricId VARCHAR(15), 
	nonSeeds VARCHAR(4500), 
	PRIMARY KEY (patricId)
);



