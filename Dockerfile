FROM python:3.7.15
RUN apt-get update
RUN apt-get install texlive-latex-base texlive-latex-extra texlive-fonts-recommended xzdec -y
RUN rm -rf /var/lib/apt/lists/*
RUN pip install pandas numpy matplotlib scipy tifffile

WORKDIR /home