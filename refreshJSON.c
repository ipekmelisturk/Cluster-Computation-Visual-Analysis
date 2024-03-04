#include <stdio.h>

int main() {
    // File path
    char file_path[] = "output.json";

    // Open file in write mode, which creates an empty file, effectively deleting its content
    FILE *file = fopen(file_path, "w");

    if (file == NULL) {
        printf("Failed to open file for writing\n");
        return 1;
    }

    // Close the file
    fclose(file);

    printf("File content deleted successfully\n");

    return 0;
}