# sndg-web
SNDG Web Apps


ALTER TABLE `sndg5`.`bioentry` 
CHANGE COLUMN `index_updated` `index_updated` TINYINT(1) NOT NULL DEFAULT 0 ;


/opt/solr-7.3.1/bin/solr start
./manage.py shell_plus --ipython --print-sql


pip install  git+https://github.com/ezequieljsosa/elsapy.git@complete_view_fields


HAYSTACK_ID_FIELD="item.id" HAYSTACK_DJANGO_CT_FIELD='metadata.django_ct'  HAYSTACK_DJANGO_ID_FIELD='item.id' HAYSTACK_DOCUMENT_FIELD='metadata.search'  ./manage.py rebuild_index -u oai


