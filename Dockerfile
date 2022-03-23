FROM ubuntu AS builder
# suppress input
ARG DEBIAN_FRONTEND=noninteractive

# essential, cmake, git, boost, GMP for Bignum library
RUN apt-get update && apt-get -y install build-essential cmake git libboost-all-dev libgmp-dev

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