#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>

int main(int argc, char **argv)
{
    int n;
    DIR *folder;
    struct dirent *dir;
    char *image_folder = "images/";
    char path[1024];
    char command[100];

    if ((folder = opendir(argv[1])) == NULL)
    {
        printf("Could not open directory: %s\n", argv[1]);
        return 1;
    }

    printf("Enter the number of ranks: ");
    scanf("%d",&n);

    sprintf(command, "./refreshJSON");
    system(command);

    while ((dir = readdir(folder)) != NULL)
    {
        if (strcmp(dir->d_name, ".") != 0 && strcmp(dir->d_name, "..") != 0)
        {
            // sprintf(command, "srun -n %d --mpi=pmix ./sobel_coarse %s%s", n, image_folder, dir->d_name);
            sprintf(command, "mpiexec -n %d ./sobel_coarse %s%s", n, image_folder, dir->d_name);

            system(command);
        }
    }
    sprintf(command, "./addBrackets");
    system(command);
    
    return 0;
}