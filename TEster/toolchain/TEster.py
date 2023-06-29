
import click
import sys
from TEster.init.sequence_generator import sequence_generator, EmptyInputFileException
from TEster.init.run_finder import run_finder
from TEster.parametrization.parameter_tester import run_analysis
from TEster.utils.tester_utils import reset_config, set_config_to_final

tool_used = None


@click.command()
@click.argument("input_file", required=True, type=click.Path(exists=True))
@click.option("element_percentage", "-p", default=70, help="Percentage of TE content")
@click.option("analysis_out_dir", "-d", default="analysis/", help='Output directory')
@click.option("sensitivity", "-s", default=200, help="Number of iterations for TEster to run")
@click.option("sequence_database", "-i", type=click.Path(exists=True), help='Reference element database')
@click.option("te_recognition_tool", "-t", default="ltr_finder", help='Specifies which recognition tool is to be parametrised, options: \"ltr_finder\", \"ltr_harvest\"')
def main(input_file, element_percentage, analysis_out_dir, sensitivity, sequence_database):

    tool_used = te_recognition_tool

    # generate from default or given database
    if not sequence_database:
        print("Database not provided, creating artificial database")
        reset_config()
        if tool_used == "ltr_finder":
            sequence_database = run_finder(input_file)
        else:
            sequence_database = run_harvest(input_file)
    try:
        element, generated_file = sequence_generator(input_file, sequence_database, element_percentage)

    # if the query sequence is empty/invalid
    except EmptyInputFileException as ex:
        print("Error: Invalid reference database file given {}".format(ex.message))
        print("Provide a valid database or run the program without a reference database")
        sys.exit(1)

    # runs analysis to detect a good configuration
    good_values, bad_values = run_analysis(generated_file, sensitivity, element, analysis_out_dir)

    # chooses best configuration
    if len(good_values["o"]) != 0:
        set_config_to_final(good_values, input_file, analysis_out_dir)
    else:
        set_config_to_final(bad_values, input_file, analysis_out_dir)



if __name__ == "__main__":
    main()
