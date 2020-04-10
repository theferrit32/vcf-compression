import os

class Config:
    def __init__(self, tabix_command_path, vcfc_command_directory, bgzip_command_path):
        self.tabix_command_path = tabix_command_path
        self.vcfc_command_directory = vcfc_command_directory
        self.bgzip_command_path = bgzip_command_path

    def get_vcfc_release_cmd(self):
        return os.path.join(self.vcfc_command_directory, 'main_release')

    def get_vcfc_timing_cmd(self):
        return os.path.join(self.vcfc_command_directory, 'main_timing')

    def get_vcfc_debug_cmd(self):
        return os.path.join(self.vcfc_command_directory, 'main_debug')
