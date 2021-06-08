# ${1} <runtimecontainername> example: RC
# ${2} <imagename> example: jonovox/nextflowcentos:latest

docker run \
  -it \
  --rm \
  --name ${1} \
  --mount type=bind,source=${PWD},target=/workflow \
  ${2} /bin/bash