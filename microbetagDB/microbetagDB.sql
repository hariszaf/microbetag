DROP DATABASE microbetagDB;
CREATE DATABASE microbetagDB; 
USE microbetagDB;

-- TABLE FOR THE PHENDB PREDICTIONS 
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
	thermophilicScore DECIMAL(5,4)    
);
SHOW VARIABLES LIKE "secure_file_priv";
SET GLOBAL local_infile=true;
-- You need to be root to run the following line: sudo mysql  -p --local-infile
LOAD DATA LOCAL INFILE '/var/lib/mysql-files/gtdb_phen_predictions.tsv' INTO TABLE phenDB;

-- Query to get the column names
SELECT `COLUMN_NAME` FROM `INFORMATION_SCHEMA`.`COLUMNS` WHERE `TABLE_SCHEMA`='microbetagDB' AND `TABLE_NAME`='phenDB';

SELECT * FROM phenDB WHERE gtdbId = 'GCF_004343615.1';
-- GTDBid is not present in the database 
SELECT * FROM phenDB WHERE gtdbId = 'GCA_013415845.1';
-- But the GCF exists 
SELECT * FROM phenDB WHERE gtdbId = 'GCF_013415845.1';
-- However, othere are just not there... 
SELECT * FROM phenDB WHERE gtdbId = 'GCA_900696425.1';
SELECT * FROM phenDB WHERE gtdbId = 'GCF_900696425.1';





-- TABLE FOR THE PATHWAY COMPLEMENTARITIES
CREATE TABLE pathwayCompelmenterarity(
	gtdbIdA VARCHAR(15),
	gtdbIdB VARCHAR(15),
	keggModule VARCHAR(15),
    
)



