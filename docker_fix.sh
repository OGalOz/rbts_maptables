#Docker fix
docker run -it -v /var/run/docker.sock:/run/docker.sock alpine chmod g+w /run/docker.sock

