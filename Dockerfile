FROM python:3.7.15
RUN apt-get update
RUN pip install pandas numpy matplotlib
WORKDIR /home