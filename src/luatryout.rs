use mlua::prelude::*;

pub fn main() -> LuaResult<()> {
    let lua = Lua::new();

    let code = std::fs::read_to_string("script.lua").expect("Failed to read Lua file");
    lua.load(&code).exec()?;
    Ok(())
}
