#include <direct.h>
#include "Utils/MatlabInterface.h"
#include "Utils/MatlabGMMDataExchange.h"
#include "GIF.h"

#define BUFFER_SIZE 4000
void setMatlabPath();

void main()
{
	setMatlabPath();

	int methodIndex;
	std::string objPath, vfPath;

	GIF program;
	program.run();
	system("pause");
}

void setMatlabPath()
{
	char path[BUFFER_SIZE];
	_fullpath(path, "..\\", BUFFER_SIZE);
	int index = strlen(path);
	//path[index - 4] = '\0';
	std::string matlabPath(path);
	matlabPath += "MatlabScripts\\";
	MatlabInterface::GetEngine().AddScriptPath(matlabPath.c_str());
}