from mul2py import buildtools
import argparse
from pathlib import Path

if __name__ == '__main__':
    # Parser arguments
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('file_path', help='Filepath to the result file')

    parser.add_argument('--overwrite', nargs='?', default=True,
                        help='Overwrite existing output files? Default is True')

    # Parse arguments
    arguments = parser.parse_args()
    filepath = Path(arguments.file_path)
    sim_overwrite = bool(arguments.overwrite)

    #Print arguments
    args = vars(arguments)
    print('\nGot {} arguments:'.format(len(args)))
    for arg in args:
        print('\t{} = {}'.format(arg, args[arg]))
    print('\n')

    #Check file type
    if not filepath.suffix == '.ecmat':
        raise ValueError('Can only convert "ecmat" files, got "{}"'.format(filepath.suffix))

    #Load and make signal
    signal = buildtools.make_signal(filepath)

    #Save signal as hyperspy file
    signal.save(filepath.with_suffix('.hspy'), overwrite=sim_overwrite)
