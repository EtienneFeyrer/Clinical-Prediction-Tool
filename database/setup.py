"""
Setup MariaDB using getpass for sudo password.
"""

import subprocess
import mysql.connector
import getpass
import sys


def setup_with_sudo_password():
    """Setup using sudo with password input."""

    print("=== MariaDB Setup (requires sudo) ===")

    # Get sudo password
    sudo_password = getpass.getpass("Enter your sudo password: ")

    # SQL commands
    sql_commands = """
CREATE USER IF NOT EXISTS 'annotation_user2'@'localhost' IDENTIFIED BY 'secure_pass_2025';
GRANT ALL PRIVILEGES ON *.* TO 'annotation_user2'@'localhost';
FLUSH PRIVILEGES;
CREATE DATABASE IF NOT EXISTS AnnotationCache;
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
    most_severe_consequence VARCHAR(255)
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
"""

    try:
        print("Executing setup commands...")

        # Use sudo -S to read password from stdin
        process = subprocess.Popen(
            ['sudo', '-S', 'mysql', '-u', 'root'],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        # Send sudo password + SQL commands
        stdout, stderr = process.communicate(input=f"{sudo_password}\n{sql_commands}")

        if process.returncode == 0:
            print("Setup completed successfully!")
            return True
        else:
            print(f"Setup failed: {stderr}")
            return False

    except Exception as e:
        print(f"Error during setup: {e}")
        return False


def test_connection():
    """Test the new user connection."""

    try:
        print("\nTesting connection as annotation_user...")
        connection = mysql.connector.connect(
            host='localhost',
            port=3306,
            user='annotation_user2',
            password='secure_pass_2025',
            database='AnnotationCache'
        )

        cursor = connection.cursor()
        cursor.execute("SHOW TABLES")
        tables = [table[0] for table in cursor.fetchall()]

        print(f"Connection successful!")
        print(f"   Tables created: {tables}")

        # Test table access
        cursor.execute("SELECT COUNT(*) FROM Annotation")
        ann_count = cursor.fetchone()[0]
        print(f"   Annotation table: {ann_count} rows")

        cursor.close()
        connection.close()
        return True

    except Exception as e:
        print(f"Connection test failed: {e}")
        return False


def create_simple_config():
    """Create a simple configuration for your app."""

    config_code = '''# Simple database configuration
# Copy this into your application
DATABASE_CONFIG = {
    'host': 'localhost',
    'port': 3306,
    'user': 'annotation_user2',
    'password': 'secure_pass_2025',
    'database': 'AnnotationCache'
}
# Usage example:
# import mysql.connector
# connection = mysql.connector.connect(**DATABASE_CONFIG)
'''

    try:
        with open('db_config.py', 'w') as f:
            f.write(config_code)
        print("Configuration saved to 'db_config.py'")
    except Exception as e:
        print(f"Error saving config: {e}")


if __name__ == "__main__":
    #This script will set up MariaDB for your annotation service
    #You need sudo access to create database users.

    if setup_with_sudo_password():
        if test_connection():
            create_simple_config()

            print("Your annotation service database is ready!")
            print("  host: localhost")
            print("  user: annotation_user2")
            print("  password: secure_pass_2025")
            print("  database: AnnotationCache")
        else:
            print("Setup completed but connection test failed")
    else:
        print("Setup failed - check your sudo password and MariaDB installation")
        sys.exit(1)