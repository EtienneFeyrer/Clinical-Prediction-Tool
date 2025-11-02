-- Database initialization script for MariaDB
-- This file is automatically executed when the container starts

USE AnnotationCache;

CREATE TABLE IF NOT EXISTS Annotation (
    Annotation_id INT AUTO_INCREMENT PRIMARY KEY,
    variant_id VARCHAR(255) UNIQUE NOT NULL,
    gene VARCHAR(255),
    allele_freq FLOAT,
    CADD FLOAT,
    OMIM VARCHAR(255),
    max_allele_freq FLOAT,
    ML_score FLOAT,
    most_severe_consequence VARCHAR(255),
    CLINSIG VARCHAR(500)
);
CREATE TABLE IF NOT EXISTS Transcript (
    id INT AUTO_INCREMENT PRIMARY KEY,
    variant_id VARCHAR(255) NOT NULL,
    transcript_id VARCHAR(255),
    polyphen FLOAT,
    protein_notation VARCHAR(255),
    REVEL FLOAT,
    Splice_AI FLOAT,
    Mane BOOLEAN DEFAULT FALSE,
    LOFTEE VARCHAR(255),
    impact VARCHAR(255),
    GERP FLOAT,
    cDNA_notation VARCHAR(255),
    consequences TEXT,
    FOREIGN KEY (variant_id) REFERENCES Annotation(variant_id)
);

-- Create indexes for better performance
CREATE INDEX idx_annotation_gene ON Annotation(gene);
CREATE INDEX idx_annotation_consequence ON Annotation(most_severe_consequence);
CREATE INDEX idx_transcript_variant ON Transcript(variant_id);
CREATE INDEX idx_transcript_id ON Transcript(transcript_id);
-- CREATE INDEX idx_transcript_gene ON Transcript(gene_symbol);

-- Insert sample data for testing (optional)
-- INSERT INTO Annotation (variant_id, gene, CADD, ML_score, REVEL, gnomAD_AF, max_allele_freq, OMIM, most_severe_consequence)
-- VALUES ('test_variant_1', 'BRCA1', 25.5, 0.85, 0.75, 0.001, 0.002, '113705', 'missense_variant');
