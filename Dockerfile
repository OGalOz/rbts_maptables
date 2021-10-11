FROM kbase/sdkbase2:python
MAINTAINER ogaloz@lbl.gov
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

RUN ls

RUN apt-get update && \
    apt-get install -y aptitude

RUN apt-get clean

RUN aptitude update -y && aptitude safe-upgrade -y && \
    apt-get install --yes \
    build-essential 

RUN apt-get install --yes apt-utils \
    perl5.24 \
    cpanminus 

RUN cpanm Getopt::Long \
    FindBin \
    File::stat

RUN apt-get install --yes r-base

RUN apt-get install python3

RUN pip install --upgrade pip 
RUN pip install pandas

# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
