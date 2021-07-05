The changes brought with this pull request are intended to make the process of app.py more automated. There are three new files, two python files, and a shell script.

The OrganaziedGenerator.py is intended to automatically perform every step from the lecture in order to produce the "06_bioactivity_data_3class_pIC50_pubchem_fp.csv" file and any other relevant file. Also, the file is programmed to ask the user about a protein and create a directory that will contain the resulted (CSV) files. Another functionality of the file is moving the needed file to the directory of the app.py. This is because the file was build on the assumed structure that is explained further at the end of this file. After finishing every operation, the file will call OrganaziedModleGenerator.py to perform its tasks.

The OrganaziedModleGenerator.py file is intended to generate the descriptor_list and the "pkl" files that will be used in the app.py.

The app.sh file is intended to be a file that organizes the invocation of all of the files, so the user will only need to work with one file. This file will invoke OrganaziedGenerator.py then will activate the condo environment "bioactivity", then will invoke the app.py file with Streamlit package. 

NOTE: Every graph was excluded from this automation.

NOTE: The assumed structure is as follow: The three files and every file that was used in part 1 to part 6, that were essential to complete the tasks of these parts, are on one directory. And inside this directory, the bioactivity-prediction-app-main directory includes everything in the original repository.