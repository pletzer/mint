FROM gcc:10 

LABEL author="bjlittle" \
      version="0.1" \
      description="MINT docker image with gcc and mesa"

#
# install system-level dependencies
#
RUN apt update -y \
  && apt upgrade -y \
  && apt install -y libgl1-mesa-dev 
