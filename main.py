from fastapi import FastAPI
from fastapi.responses import JSONResponse

from pydantic import BaseModel
import os

limbstager_exe = "/Users/lauavino/Documents/Code/staging-system-api/limb-staging/src/limbstager"

app = FastAPI()


def grep(filename, search_word):
    with open(filename, "r") as f:
        for line in f:
            if search_word in line:
                return line.strip().split()
    return []


class PointsData(BaseModel):
    points: list
    header: str


@app.get("/")
async def root():
    return {"message": "Welcome !!"}


@app.post("/stage/")
async def stage(pointsData: PointsData):
    print(pointsData)

    # Create a temporary file with the measurments
    measure_file = ".tmp_input.txt"
    output_file = ".tmp_output.txt"
    with open(measure_file, "w") as f:
        f.write(pointsData.header)
        for p in pointsData.points:
            f.write(f"MEASURED {p[0]} {p[1]}\n")

    # Check that the executable exits
    if not os.path.isfile(limbstager_exe):
        return JSONResponse(status_code=418)

    # now stage: a .tmp_out.txt file is created
    errnr = os.system(f"{limbstager_exe} {measure_file} > {output_file}")

    if errnr:
        print(f"limbstager executable {limbstager_exe} returned error:", errnr)
        return {"message": "Oh oh"}
    else:
        print(f"Response of the limbstage: {errnr}")

    # print("Reasing the grep!")
    result = grep(output_file, "RESULT")
    # print(result, "aaa")

    if not len(result):
        # print("Error - Could not stage the limb, RESULT tag is missing")
        return JSONResponse(
            status_code=418,
            content={
                "message":
                "Error - Could not stage the limb, RESULT tag is missing"
            },
        )
    stage = result[1]
    print(stage)
    # txt.text(f"Limb staged as {stage}")
    # plt.at(0).render()

    return JSONResponse(
        status_code=200,
        content={
            "message": "Stage succesfully computed",
            "stage": stage
        },
    )
