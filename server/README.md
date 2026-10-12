# microbetag back-end

## Basic architecture and firing up

microbetag is hosted in KU Leuven server and is Docker oriented. 
More specifically, it consists of 3 Docker images:
- a database (`services/db/`)
- an ngninx (`services/nginx`)
- the actual app (`services/web`)
 
First, each of those images needs to be built, and their settings filled in: copy `.env.example` to `.env.dev`,
and `services/db/.env.example` to `services/db/.env`, and set the values (they are not tracked by git).
Then `docker compose up`, run from this folder, starts the web server.

### Running it as a service (systemd)

To keep the server running after you log out, start it at boot and restart it if it fails, let systemd run it.
[`deploy/microbetag.service`](./deploy/microbetag.service) is a unit file for this; set its `WorkingDirectory`
to this folder (and the path of `docker`, from `which docker`, if it is not `/usr/bin/docker`), then:

```bash
sudo cp deploy/microbetag.service /etc/systemd/system/
sudo systemctl daemon-reload
sudo systemctl enable --now microbetag
```

`enable` starts it at every boot, `--now` starts it right away.

| Task | Command |
|---|---|
| Is it running? | `systemctl status microbetag` |
| Follow the logs | `journalctl -u microbetag -f` |
| Restart, e.g. after rebuilding an image | `sudo systemctl restart microbetag` |
| Stop | `sudo systemctl stop microbetag` |
| Do not start it at boot | `sudo systemctl disable microbetag` |

systemd restarts the whole stack if `docker compose up` fails. To restart a single crashed container without
the others, add `restart: unless-stopped` to its service in `docker-compose.yml`.

The app is now up and running!
You may check this by just typing to a browser: [msysbio.gbiomed.kuleuven.be/](msysbio.gbiomed.kuleuven.be/).


## How to edit the app

As mentioned above, the app consists of 3 Docker images that all 3 of them are used simultaneously.
If one is to edit something on the app, they need to move to the corresponding 
image, edit the code and then **rebuild** the Docker image. 
For example, if you are about to edit a route of the `app.py` script, then you will have to rebuild the `microbetagapp-web` image:

```bash
docker build -f Dockerfile.prod -t microbetagapp-web:latest .
```

To make sure you can always come back to a working version of the app, a good pracice is to first make a backup **tag** of your working version
running

```bash
docker tag <docker_image_id> microbetagapp-web-working-version
```

and only then rebuild the Docker image. 


Once you have rebuilt the image, restart the service:

```bash
sudo systemctl restart microbetag
```

This re-creates the `web` container, and after a few seconds the app is up and running again with your new code.



## Notes on the server's environment


## Env

I made a new file under  `/etc/docker/` called `daemon.json`:

```bash
u0156635@gbw-s-msysbio01:/etc/docker$ more daemon.json 
{
  "data-root": "/data/docker"

}
```


### Ports 

I changed the mysql port of the vm to 4545 so I can have the 3306 for the container


### My database

We need to build the whole database in the vm or wherever (even on a laptop) and then export it to a `dbtest_data.sql` file. 
Here is how: 

```bash
 mysqldump -u [user name] –p [password] [database_name] > [dumpfilename.sql]
```

### To build the SSL keys 

Follow the instructions [here](https://www.baeldung.com/openssl-self-signed-cert)
and remove the pass phrase as shown [here](https://help.cloud66.com/docs/security/remove-passphrase).

> `openssl rsa -in [original.key] -out [new.key]`


## Links

- [Setting up and Running a MySQL Container](https://www.baeldung.com/ops/docker-mysql-container)
- [Creating a Self-Signed Certificate With OpenSS](https://www.baeldung.com/openssl-self-signed-cert)

Challenges:
- [`argparse` in nginx](https://github.com/benoitc/gunicorn/issues/1867)



## Notes after destroying `/data`


failed to update store for object type *libnetwork.endpointCnt: Key not found in store

this fixes by running 
sudo service docker restart


## VM overall

consider of crontabs removing for example /tmp folders



## Notes Spring 2025


To set the whole thing locally for development purposes:


- Install `docker compose` 

  ```
  mkdir -p ~/.docker/cli-plugins/
  curl -SL https://github.com/docker/compose/releases/download/v2.3.3/docker-compose-linux-x86_64 -o ~/.docker/cli-plugins/docker-compose
  chmod +x ~/.docker/cli-plugins/docker-compose
  ```



  conda create -n mtg_app



🔴 Remember to first build each docker image on each own when you edit something, since docker compose will use the already built images! 


Check that locally everything is up and running on http://localhost:1337
the port that has been denoted in the main `docker-compose.yml` script, on the `nginx` container.

