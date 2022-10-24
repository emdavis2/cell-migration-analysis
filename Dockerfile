# FROM python:3.7.15
# RUN apt-get update
# RUN pip install pandas numpy matplotlib scipy tifffile

# WORKDIR /home

FROM ubuntu:bionic

RUN ln -snf /usr/share/zoneinfo/Etc/UTC /etc/localtime \
    && echo "Etc/UTC" > /etc/timezone \
    && apt-get update \
    && apt-get upgrade -y \
    && apt-get install texlive-latex-base texlive-latex-extra texlive-fonts-recommended xzdec -y \
    && rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y software-properties-common gcc && \
    add-apt-repository -y ppa:deadsnakes/ppa

RUN apt-get update && apt-get install -y python3.7 python3-distutils python3-pip python3-apt

#RUN pip3 install Cython pandas numpy matplotlib scipy tifffile

WORKDIR /home