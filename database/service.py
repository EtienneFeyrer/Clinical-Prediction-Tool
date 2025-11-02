import os
import mysql.connector


class DatabaseService:
    """Database operations for annotations."""

    @staticmethod
    def get_db_config():
        """Get database configuration from environment variables."""
        return {
            'host': os.getenv('DB_HOST', 'localhost'),
            'port': int(os.getenv('DB_PORT', '3306')),
            'user': os.getenv('DB_USER', 'annotation_user'),
            'password': os.getenv('DB_PASSWORD', 'secure_pass_2024'),
            'database': os.getenv('DB_NAME', 'AnnotationCache')
        }

    @staticmethod
    def variant_exists(variant_key: str) -> bool:
        """Check if variant already exists in database."""
        try:
            connection = mysql.connector.connect(**DatabaseService.get_db_config())
            cursor = connection.cursor()
            cursor.execute("SELECT COUNT(*) FROM Annotation WHERE variant_id = %s", (variant_key,))
            result = cursor.fetchone()
            cursor.close()
            connection.close()
            return result[0] > 0

        except mysql.connector.Error as e:
            print(f"Database error in variant_exists: {e}")
            return False
        except Exception as e:
            print(f"Unexpected error in variant_exists: {e}")
            return False

    @staticmethod
    def get_variant_annotation(variant_key: str) -> dict:
        """Get annotation data for variant."""
        try:
            connection = mysql.connector.connect(**DatabaseService.get_db_config())
            cursor = connection.cursor()

            # Get main annotation data
            main_query = """
                         SELECT gene, \
                                CADD, \
                                ML_score, \
                                most_severe_consequence,
                                allele_freq, \
                                max_allele_freq, \
                                OMIM,
                                CLINSIG
                         FROM Annotation \
                         WHERE variant_id = %s \
                         """
            cursor.execute(main_query, (variant_key,))
            main_row = cursor.fetchone()

            if not main_row:
                cursor.close()
                connection.close()
                return None

            # Get transcript data
            transcript_query = """
                               SELECT transcript_id, \
                                      polyphen, \
                                      protein_notation, \
                                      REVEL, \
                                      Splice_AI,
                                      Mane, \
                                      LOFTEE, \
                                      impact, \
                                      GERP, \
                                      cDNA_notation, \
                                      consequences
                               FROM Transcript \
                               WHERE variant_id = %s \
                               """
            cursor.execute(transcript_query, (variant_key,))
            transcript_rows = cursor.fetchall()

            cursor.close()
            connection.close()

            # Build response
            transcripts = []
            for row in transcript_rows:
                transcript = {
                    'transcript_id': row[0],
                    'polyphen': row[1],
                    'protein_notation': row[2],
                    'REVEL': row[3],
                    'Splice_AI': row[4],
                    'Mane': bool(row[5]) if row[5] is not None else False,
                    'LOFTEE': row[6],
                    'impact': row[7],
                    'GERP': row[8],
                    'cDNA_notation': row[9],
                    'consequences': row[10]
                }
                transcripts.append(transcript)

            return {
                'gene': main_row[0],
                'CADD': main_row[1],
                'pathogenicity_score': main_row[2],  # ML_score
                'Most Severe Consequence': main_row[3],
                'gnomAD AF': main_row[4],  # allele_freq
                'max_allele_freq': main_row[5],
                'OMIM': main_row[6],
                'CLINSIG': main_row[7],
                'transcript_consequences': transcripts
            }

        except mysql.connector.Error as e:
            print(f"Database error in get_variant_annotation: {e}")
            return None
        except Exception as e:
            print(f"Unexpected error in get_variant_annotation: {e}")
            return None

    @staticmethod
    def bulk_insert_annotations(annotations_data, transcripts_data):
        """
        Bulk insert annotations and transcripts into database.

        Args:
            annotations_data: List of tuples for Annotation table
            transcripts_data: List of tuples for Transcript table

        Returns:
            True if successful, False otherwise
        """
        connection = None
        try:
            connection = mysql.connector.connect(**DatabaseService.get_db_config())
            cursor = connection.cursor()

            # Bulk insert annotations
            cursor.executemany("""
                               INSERT INTO Annotation (variant_id, gene, CADD, ML_score, most_severe_consequence,
                                                       allele_freq, max_allele_freq, OMIM, CLINSIG)
                               VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s) ON DUPLICATE KEY
                               UPDATE CADD =
                               VALUES (CADD)
                               """, annotations_data)

            # Bulk insert transcripts
            if transcripts_data:
                cursor.executemany("""
                                   INSERT INTO Transcript (variant_id, transcript_id, polyphen, protein_notation,
                                                           REVEL, Splice_AI, Mane, LOFTEE, impact, GERP,
                                                           cDNA_notation, consequences)
                                   VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
                                   """, transcripts_data)

            connection.commit()
            cursor.close()
            connection.close()
            print(f"BULK STORED {len(annotations_data)} variants successfully")
            return True

        except Exception as e:
            print(f"Bulk storage error: {e}")
            if connection:
                connection.rollback()
                connection.close()
            return False
