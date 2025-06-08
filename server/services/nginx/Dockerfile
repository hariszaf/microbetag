FROM nginx:1.25

RUN apt-get update && apt-get install -y apt-transport-https
RUN apt-get install -y vim

RUN rm /etc/nginx/conf.d/default.conf
COPY nginx.conf /etc/nginx/conf.d
RUN sed -E -i '17i \ \ \ \ client_max_body_size 20M;' /etc/nginx/nginx.conf

# The proxy_read_timeout in Nginx specifies the maximum time Nginx will wait for a response from the backend server (Gunicorn). 
RUN sed -E -i '30i \ \ \ \ proxy_read_timeout 600;' /etc/nginx/nginx.conf
RUN sed -E -i '31i \ \ \ \ proxy_connect_timeout 600;' /etc/nginx/nginx.conf
RUN sed -E -i '32i \ \ \ \ proxy_send_timeout 600;' /etc/nginx/nginx.conf

RUN sed -E -i '34i \ \ \ \ fastcgi_read_timeout 600;' /etc/nginx/nginx.conf

RUN sed -E -i '36i \ \ \ \ proxy_set_header Connection "";' /etc/nginx/nginx.conf
