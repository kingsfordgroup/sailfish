execute_process(COMMAND tar xzvf sample_data.tgz
                WORKING_DIRECTORY ${TOPLEVEL_DIR}
                RESULT_VARIABLE TAR_RESULT
               )
execute_process(COMMAND ${TOPLEVEL_DIR}/bin/sailfish index -t transcripts.fasta -k 20 -o sample_index
                WORKING_DIRECTORY ${TOPLEVEL_DIR}/sample_data
                RESULT_VARIABLE INDEX_RESULT
                )
execute_process(COMMAND ${TOPLEVEL_DIR}/bin/sailfish quant -i sample_index -r reads_1.fastq reads_2.fastq -o sample_quant
	            WORKING_DIRECTORY ${TOPLEVEL_DIR}/sample_data
	            RESULT_VARIABLE QUANT_RESULT
                )
if (EXISTS ${TOPLEVEL_DIR}/sample_data/sample_quant/quant.sf)
	message("Sailfish ran successfully")
else()
	message("Sailfish failed to produce output")
endif()