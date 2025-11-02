# Development Setup with Docker Compose

This document explains how to set up the annotation server for local development using Docker Compose.

## Prerequisites

- Docker
- Docker Compose
- install docker with: 
```bash 
sudo apt install docker.io docker-compose-plugin
```


## Quick Start

1. **Start the services:**
   ```bash
   docker compose up -d
   ```

2. **Check service status:**
   ```bash
   docker compose ps
   ```

3. **View logs:**
   ```bash
   # All services
   docker compose logs -f
   
   # Specific service
   docker compose logs -f api
   docker compose logs -f mariadb
   ```

4. **Test the API:**
   ```bash
   curl http://localhost:5001/health
   curl http://localhost:5001/statistics
   ```

5. **Stop the services:**
   ```bash
   docker compose down
   ```

## Database Access

### Connect to MariaDB
```bash
# Using docker exec
docker-compose exec mariadb mysql -u annotation_user -psecure_pass_2024 AnnotationCache

# Using mysql client from host (if installed)
mysql -h 127.0.0.1 -P 3307 -u annotation_user -psecure_pass_2024 AnnotationCache
```

### Database Management
```bash
# Reset database (removes all data)
docker compose down -v
docker compose up -d

# View database tables
docker compose exec mariadb mysql -u annotation_user -psecure_pass_2024 AnnotationCache -e "SHOW TABLES;"

#View Annotation Table:
docker compose exec mariadb mysql -u annotation_user -psecure_pass_2024 AnnotationCache -e "SELECT * FROM Annotation;"


# Truncate tables for clean testing
docker compose exec mariadb mysql -u annotation_user -psecure_pass_2024 AnnotationCache -e "SET FOREIGN_KEY_CHECKS = 0; TRUNCATE TABLE Transcript; TRUNCATE TABLE Annotation; SET FOREIGN_KEY_CHECKS = 1;"
```

## Development Workflow

1. **Make code changes** - Files are mounted as volumes, so changes are reflected immediately
2. **Restart API service** if needed:
   ```bash
   docker compose restart api
   ```
3. **View real-time logs:**
   ```bash
   docker compose logs -f api
   ```

## Configuration

### Environment Variables
The services use the following environment variables:

- `DB_HOST`: MariaDB hostname (mariadb)
- `DB_USER`: Database user (annotation_user)
- `DB_PASSWORD`: Database password (secure_pass_2024)
- `DB_NAME`: Database name (AnnotationCache)

### Ports
- **API Server**: http://localhost:5001
- **MariaDB**: localhost:3306

## Troubleshooting

### API Server won't start
```bash
# Check if MariaDB is healthy
docker-compose exec mariadb mysqladmin ping --silent

# Check API logs
docker-compose logs api
```

### Database connection issues
```bash
# Verify MariaDB is running
docker-compose ps mariadb

# Test database connection
docker-compose exec mariadb mysql -u annotation_user -psecure_pass_2024 AnnotationCache -e "SELECT 1;"
```

### Reset everything
```bash
docker compose down -v --remove-orphans
docker compose up -d
```