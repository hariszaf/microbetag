DROP DATABASE microbetagDB;
CREATE DATABASE microbetagDB;
USE microbetagDB;

DROP TABLE Pathway_complementarity; -- you first need to drop the children tables 
DROP TABLE Identifiers;


CREATE TABLE IF NOT EXISTS Identifiers(
--     Include: KEGG terms, GENOME ids, NCBI Tax Ids, 
    term VARCHAR(15),
    term_description VARCHAR(50),
    ontology VARCHAR(10),
    url VARCHAR(100),
    PRIMARY KEY(term)
);

INSERT INTO Identifiers VALUES ("mboi", "kegg genome", "KEGG", "https://www.kegg.jp/entry/K00001");
INSERT INTO Identifiers VALUES ("GCA_002006445.1", "GTDB genome", "RefSeq", "end.com");
INSERT INTO Identifiers VALUES ("476",  "NCBI Tax Id", "NCBI", "476end.com");
INSERT INTO Identifiers VALUES ("1520", "NCBI Tax Id", "NCBI", "1520end.com");

SELECT * FROM Identifiers;


--     kegg_module_url VARCHAR(100),    this will be provdided through the identifiers table
CREATE TABLE Pathway_complementarity(
    beneficiarys_ncbiId VARCHAR(8),
    beneficiarys_genome VARCHAR(15),
    donors_ncbiId VARCHAR(8),
    donors_genome VARCHAR(15),
    kegg_module VARCHAR(6),
    complement VARCHAR(50),
    pathway VARCHAR(150),
    complement_percentage FLOAT(3),
	FOREIGN KEY (beneficiarys_ncbiId) REFERENCES Identifiers (term), -- ON UPDATE CASCADE ON DELETE CASCADE,
	FOREIGN KEY (beneficiarys_genome) REFERENCES Identifiers (term), -- ON UPDATE CASCADE ON DELETE CASCADE,
	FOREIGN KEY (donors_ncbiId) REFERENCES Identifiers (term), -- ON UPDATE CASCADE ON DELETE CASCADE,
	FOREIGN KEY (donors_genome) REFERENCES Identifiers (term) -- ON UPDATE CASCADE ON DELETE CASCADE    
);

INSERT INTO Pathway_complementarity (beneficiarys_ncbiId, beneficiarys_genome, donors_ncbiId, donors_genome, kegg_module, complement, pathway, complement_percentage ) 
VALUES ("476", "mboi", "1520", "GCA_002006445.1", "M00096", "K01823", "K01662K00099|K00991|K00919|K01770|K03526|K03527|K01823", 1/8);

SELECT * FROM Pathway_complementarity;
SELECT * FROM Identifiers WHERE term =  (SELECT donors_ncbiId FROM Pathway_complementarity WHERE donors_ncbiId = 1520);

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
);

SHOW VARIABLES LIKE "secure_file_priv";
SET GLOBAL local_infile=true;
-- You need to be root to run the following line: sudo mysql  -p --local-infile
-- LOAD DATA LOCAL INFILE '/var/lib/mysql-files/gtdb_phen_predictions.tsv' INTO TABLE phenDB;



SELECT * FROM phenDB WHERE gtdbId = "GCF_910593785.1";

SELECT  rAcetoin  FROM phenDB WHERE gtdbId = "GCF_910593785.1";


-- --  ALTER TABLE phenDB ENGINE=InnoDB/MyISAM;
-- SHOW VARIABLES LIKE "secure_file_priv";
-- SET GLOBAL local_infile=true;

