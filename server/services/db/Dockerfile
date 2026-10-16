# The microbetag database: MySQL, with the schema of every microbetag database. At its first start (an empty data
# folder), MySQL runs schema/microbetagDB.sql, which creates the tables and loads the files of microbetag-export-db
# from /var/lib/mysql-files (mounted, see docker-compose.yml).
#
# The credentials (MYSQL_DATABASE, MYSQL_USER, MYSQL_PASSWORD, MYSQL_ROOT_PASSWORD) are given when the container
# starts (env_file of docker-compose.yml), not built into the image.

FROM mysql:8.4

COPY schema/microbetagDB.sql /docker-entrypoint-initdb.d/

EXPOSE 3306
