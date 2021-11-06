docker run --name fasterq -ti -d -v $PWD:/data staphb/sratoolkit
docker exec -ti fasterq ./download.sh
docker stop fasterq
docker rm fasterq
