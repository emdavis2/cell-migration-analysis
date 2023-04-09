FROM python:3.9.16
RUN apt-get update
#RUN apt-get install texlive-latex-base texlive-latex-extra texlive-fonts-recommended xzdec -y
#RUN rm -rf /var/lib/apt/lists/*
RUN pip install pandas==1.4.4 numpy matplotlib scipy tifffile seaborn

WORKDIR /home