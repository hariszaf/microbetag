
FROM mysql:latest

RUN microdnf install -y vim
RUN microdnf install --assumeyes python3-pip
RUN pip3.9 install mysql-connector-python

RUN chown -R mysql:root /var/lib/mysql/ 

# The credentials (MYSQL_DATABASE, MYSQL_USER, MYSQL_PASSWORD, MYSQL_ROOT_PASSWORD) are given when the container
# starts (env_file of docker-compose.yml), not built into the image

COPY requirements.txt .
RUN pip install -r requirements.txt

# COPY init.sql /docker-entrypoint-initdb.d/
COPY schema/local_init.sql /docker-entrypoint-initdb.d/
COPY scripts/ /docker-entrypoint-initdb.d/scripts/

# NOTE: We do not COPY the data themselves, as we have a symbolic link to them pointing to /var/lib/mysql-files

EXPOSE 3306
