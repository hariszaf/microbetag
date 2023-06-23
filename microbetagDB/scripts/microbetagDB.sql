DROP DATABASE microbetagDB;
CREATE DATABASE microbetagDB; 
USE microbetagDB;

-- ///////////////////////////////////////
-- TABLE RELATADE TO PHEND PREDICTIONS
-- ///////////////////////////////////////

DROP TABLE phenDB; 
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
    isovalericAcidScire DECIMAL(5,4), 
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
    PRIMARY KEY (gtdbId)
)  
ENGINE=InnoDB 
ROW_FORMAT=COMPRESSED;

-- Disable keys and indexes
ALTER TABLE phenDB DISABLE KEYS;

-- You need to be root to run the following line: sudo mysql  -p --local-infile
LOAD DATA INFILE '/var/lib/mysql-files/gtdb_phen_predictions.tsv'
INTO TABLE phenDB
FIELDS TERMINATED BY '\t'
LINES TERMINATED BY '\n'
IGNORE 1 LINES;

-- Enable again keys
ALTER TABLE phenDB ENABLE KEYS;

-- Build indexes for the gtdb ids 
CREATE INDEX idx_column ON phenDB (gtdbId);
ANALYZE TABLE phenDB;

-- ///////////////////////////////////////
-- PATHWAY COMPLEMENTARITY RELATED TABLES 
-- ///////////////////////////////////////

-- Query to get the column names
SELECT `COLUMN_NAME` FROM `INFORMATION_SCHEMA`.`COLUMNS` WHERE `TABLE_SCHEMA`='microbetagDB' AND `TABLE_NAME`='phenDB';



-- Table for GenomeId 2 NCBI TaxId 
CREATE TABLE genome2taxNcbiId(
	genomeId VARCHAR(15),
    ncbiTaxId VARCHAR(10),
    PRIMARY KEY (genomeId)
);

-- Table for unique complementarity cases
CREATE TABLE uniqueComplements(
	complementId INT AUTO_INCREMENT,
    KoModuleId VARCHAR(7),
    complement VARCHAR(50), 
    pathway VARCHAR(250),
    colorMapUrl VARCHAR(1000),  -- CONSIDER REMOVING THIS
    PRIMARY KEY (complementId)
);


-- TABLE FOR THE PATHWAY COMPLEMENTARITIES
DROP TABLE pathwayComplementarity;
CREATE TABLE pathwayComplementarity(
	beneficiaryGenome VARCHAR(15),
	donorGenome VARCHAR(15),
	complmentId INT(4),
    FOREIGN KEY (beneficiaryGenome) REFERENCES genome2taxNcbiId(genomeId),
    FOREIGN KEY (donorGenome) REFERENCES genome2taxNcbiId(genomeId),
    FOREIGN KEY (complmentId) REFERENCES uniqueComplements(complementId)
);
ALTER TABLE pathwayComplementarity ADD PRIMARY KEY(beneficiaryGenome, donorGenome);


INSERT INTO genome2taxNcbiId VALUES("GCA_004365965.1", "364297");
INSERT INTO genome2taxNcbiId VALUES("tlo", "2891210");
INSERT INTO uniqueComplements VALUES (1, "M00134", "K01476", "K01476;K01581", "urllll");
INSERT INTO pathwayComplementarity VALUES ("GCA_004365965.1", "tlo", 1);

-- Get all the complementarity entries of a specific gneome
SELECT  *  FROM pathwayComplementarity WHERE beneficiaryGenome = "GCA_004365965.1";







-- Get all pathway complementarities for all genomes of the beneficary
SELECT beneficiaryGenome FROM pathwayComplementarity
INNER JOIN genome2taxNcbiId 
ON pathwayComplementarity.beneficiaryGenome = genome2taxNcbiId.genomeId
WHERE genome2taxNcbiId.ncbiTaxId = "364297" ;


-- Get all pathway complementarities for all genomes of the donor
SELECT donorGenome FROM pathwayComplementarity
INNER JOIN genome2taxNcbiId 
ON pathwayComplementarity.donorGenome = genome2taxNcbiId.genomeId
WHERE genome2taxNcbiId.ncbiTaxId = "2891210" ;


-- Get ALL the complementarities between 2 genomes - OK
SELECT uniqueComplements.KoModuleId, uniqueComplements.complement, uniqueComplements.pathway
FROM uniqueComplements
INNER JOIN  pathwayComplementarity
ON pathwayComplementarity.complmentId = uniqueComplements.complementId
WHERE  pathwayComplementarity.beneficiaryGenome = "GCA_004365965.1" AND pathwayComplementarity.donorGenome = "tlo";





-- Get all complements ids for a specific NCBI Taxonomy Id // we do not care about which genome ids are used !
-- WE STILL NEED TO HAVE A WHERE FOR THE DONOR'S NCBI ID!
SELECT 
	complmentId FROM pathwayComplementarity
RIGHT JOIN 
	genome2taxNcbiId 
ON 
	pathwayComplementarity.beneficiaryGenome = genome2taxNcbiId.genomeId 
WHERE 
	genome2taxNcbiId.ncbiTaxId = "364297" ;
