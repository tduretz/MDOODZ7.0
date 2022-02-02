import subprocess
import h5py

mdoodz_source_dir = 'SOURCE'


class MdoodzModel:
    model_name = str
    argv = [str]

    def __init__(self, model_name: str):
        self.model_name = model_name

    def compile(self):
        stream = subprocess.Popen(
            ['make', 'clean', 'all', f'MODEL={self.model_name}', 'OPT=no', 'MKL=no'], cwd=mdoodz_source_dir)
        print(stream.stdout)
        stream.wait()

    def set_params(self, argv: [str]):
        self.argv = argv
        # TODO set arguments

    def run(self):
        stream = subprocess.Popen([
            f'./Doodzi_{self.model_name}',
            f'{self.model_name}.txt',
        ], cwd=mdoodz_source_dir)
        print(stream.stdout)
        stream.wait()

    def get_results(self):
        file_name = 'Output00001.gzip.h5'
        return h5py.File(f'{mdoodz_source_dir}/{file_name}', 'r')
