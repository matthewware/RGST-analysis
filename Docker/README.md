# RGST Docker

This Dockerfile creates a docker image with the the suite of software needed to reproduce the analysis for the randomized GST experiments. These instructions assume you have a working install of Docker.

To build, simply run the following from the top level of this repo.

```bash
docker build -t rgst_env -f Docker/Dockerfile .
```

To run things there's a convenience scripts:

```bash
./start_env.sh
```

This will print a URL that you can paste into your browser of choice.

In the container is a `notebook/` folder with the jupyter analysis notebooks.  
Add any files or figures you want to save to the `things-to-save` folder and
run the `copy_data.sh` script to get them back to your host machine.  To shut
things done, close all the browser tabs and run `stop_env.sh`.

### Trouble shooting

If you have trouble with any of the scripts, you can always run things manually:

```bash
$ docker ps
CONTAINER ID        IMAGE               COMMAND                  CREATED             STATUS              PORTS                    NAMES
1b2970bb07a4        rgst_env:latest     "tini -- start-noteb…"   8 seconds ago       Up 7 seconds        0.0.0.0:8888->8888/tcp   exciting_liskov
$ docker port 1b2970bb
8888/tcp -> 0.0.0.0:32781
$ docker logs --tail 3 1b2970bb
Copy/paste this URL into your browser when you connect for the first time,
   to login with a token:
       http://localhost:8888/?token=9d1543e0f50b0097f703d0abeae4f8ed71b692caf542c654  
```

and point your browser to `http://localhost:32781/?token=9d1543e0f50b0097f703d0abeae4f8ed71b692caf542c654`.

You should now have a working kernel for python2, python3 and julia v0.6.3
