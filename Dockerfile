FROM ubuntu AS builder
# suppress input
ARG DEBIAN_FRONTEND=noninteractive

# essential, cmake, git, boost, GMP for Bignum library
RUN \
apt-get update -y && \
apt-get upgrade -y && \
apt-get -y install build-essential git libboost-all-dev libgmp-dev wget && \
apt-get autoremove -y && \
apt-get clean -y

ARG CMAKE_INSTALLER=cmake-3.22.3-linux-x86_64.sh
RUN \
wget https://github.com/Kitware/CMake/releases/download/v3.22.3/${CMAKE_INSTALLER} && \
chmod +x ${CMAKE_INSTALLER}

RUN ./${CMAKE_INSTALLER} --prefix=/usr/local --exclude-subdir

# unittest-cpp
RUN git clone https://github.com/unittest-cpp/unittest-cpp.git
WORKDIR /unittest-cpp/builds
RUN cmake ../
RUN cmake --build ./
RUN cmake --build ./ --target install

# libalgebra_tests build
WORKDIR /libalgebra
COPY . .
WORKDIR /libalgebra/build
ENTRYPOINT cmake -DCMAKE_BUILD_TYPE=Release -DLIBALGEBRA_TESTING=ON .. && cmake --build .

# Run the tests
FROM builder AS tester
ENTRYPOINT cmake -DCMAKE_BUILD_TYPE=Release -DLIBALGEBRA_TESTING=ON .. && cmake --build . && ctest

# docker volume create libalgebra_volume
# docker build -t libalgebra . --target builder
# docker run -v libalgebra_volume:/libalgebra/build -it libalgebra
#
# OR
#
# docker volume create libalgebra_tests_volume
# docker build -t libalgebra_tests . --target tester
# docker run -v libalgebra_tests_volume:/libalgebra/build -it libalgebra_tests
