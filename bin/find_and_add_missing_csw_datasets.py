import csv
import sys
import os
from geophys_utils import CSWUtils

csv_from_oracle = sys.argv[1]
csv_from_csw_find = sys.argv[2]
# Python code t get difference of two lists
# Using set()
output_csv_path = sys.argv[3]


def find_missing_datasets():
    """
    Iterates through two csv files and finds matches the datasets based on the "file_path" column in the csv derived
     from the oracle database and "url" from the csw_find csv.
    :return: A list of the datasets urls that are in the oracle csv, but not the csw_find csv
    e.g. ['/g/data/uc0/rr2_dev/axi547/GADDS/non_awags_line/GSV3023MAG.nc', '/g/data/uc0/rr2_dev/axi547/GADDS/non_awags_line/P652DEM.nc']
    """
    with open(csv_from_oracle, "r") as datasets_from_oracle_csv:
        csv_reader_datasets = csv.DictReader(datasets_from_oracle_csv, delimiter=',')
        datasets_from_orcale_file_paths = []
        for lines in csv_reader_datasets:
            datasets_from_orcale_file_paths.append(lines['file_path'])

    with open(csv_from_csw_find, "r") as csw_find_csv:
        csv_reader_csw = csv.DictReader(csw_find_csv, delimiter=',')
        csw_find_result_urls = []
        for lines in csv_reader_csw:
            csw_find_result_urls.append(lines['url'])

    datasets_csw_failed_to_find = list(set(datasets_from_orcale_file_paths) - set(csw_find_result_urls))
    return datasets_csw_failed_to_find


def get_dataset_title_of_missing_datasets(missing_datasets_list):
    """
    Finds the dataset title in the csv from the list of urls
    :param missing_datasets_list:
    :return: list of titlewords
    """

    list_of_titlewords_to_search = []
    with open(csv_from_oracle, "r") as datasets_from_oracle_csv:
        csv_reader_datasets = csv.DictReader(datasets_from_oracle_csv, delimiter=',')
        for lines in csv_reader_datasets:

            for url in missing_datasets_list:
                if (lines['file_path']) == url:
                    list_of_titlewords_to_search.append(lines['DATASET_TITLE'])
        return list_of_titlewords_to_search



def call_csw_query_for_each_missing_dataset_and_write_csv(list_of_titlewords):
    """
    Wrangles each list of titles words into the correct format. Words like "of" make the query fail. All words must be
    comma separated with no spaces.
    Then for each title in will do a csw_find titlewords search, and add the result to the output csv.
    :param list_of_titlewords:
    :return:
    """

    #test_list = ["Total Magnetic Intensity (TMI) grid of Huckitta, NT, 1981 survey", "Radiometric thorium equivalent grid of Maurice, SA, 1992 survey"]

    with open(output_csv_path, "w") as output_csv:
        for i in range(len(list_of_titlewords)):
        #for i in range(len(test_list)):
            """
            convert this:   Total Magnetic Intensity (TMI) grid of Huckitta, NT, 1981 survey
            to this:        Total,Magnetic,Intensity,(TMI),grid,Huckitta,NT,1981,survey
            """
            tw = test_list[i]
            print("TW")
            print(tw)
            tw = tw.replace(' of ', ' ')
            tw = tw.replace(' the ', ' ')
            tw = tw.replace(' and ', ' ')
            tw = tw.replace(' in ', ' ')
            tw = tw.replace(' for ', ' ')
            tw = tw.replace(' at ', ' ')

            tw = tw.replace(',', '')
            tw = tw.replace(' ', ',')
            print(tw)
            if i == 0: # include header
                result = os.popen('csw_find --titlewords={}, --types=dataset, --protocols=FILE:GEO --fields="url, uuid, identifiers, title" --delimiter=, --header_row'.format(tw)).read()
                output_csv.writelines(result)
            else:
                result = os.popen('csw_find --titlewords={}, --types=dataset, --protocols=FILE:GEO --fields="url, uuid, identifiers, title" --delimiter=,'.format(tw)).read()
                output_csv.writelines(result)

def main():

    missing_datasets_list = find_missing_datasets()
    print("Datasets not found by csw_find: {}".format(missing_datasets_list))

    list_of_titlewords_to_search = get_dataset_title_of_missing_datasets(missing_datasets_list)
    call_csw_query_for_each_missing_dataset_and_write_csv(list_of_titlewords_to_search)


if __name__ == '__main__':
    main()

